# ql_grid_setup.r

# set up grid for Q-learning search in lambda-T parameter space
# to be called by ql_master.r - has depends on params defined in that file

##################
# set up grid
##################

lambda_min<-0
lambda_max<-1

lambda_points<-seq(from=lambda_min, to=lambda_max, length.out=lambda_numpoints)
temp_points<-seq(from=temp_min, to=temp_max, length.out=temp_numpoints) # note: temp means temperature, not temporary

# list of length 2 containing param values
point_lists<-list(lambda_points, temp_points)
names(point_lists)<-c('lambda', 'temperature')

# matrix containing combinations of param values
point_grid<-expand.grid(point_lists)

# list of combinations of param values
point_list<-as.list(as.data.frame(t(point_grid)))

# unique string identifiers for each param combination
indexer<-apply(point_grid, 1, paste, collapse='_')

# set target temperature and key points
temp_target<-which.min(abs(temp_points-1.5))
state_init<-c(lambda_points[1], temp_points[temp_target])
state_target<-c(lambda_points[lambda_numpoints], temp_points[temp_target])
indexer_init<-paste(state_init, collapse='_')
indexer_target<-paste(state_target, collapse='_')

################
# define intermediate distributions
################

# mu update rule for normals:
mu_update<-function(lambda, T) {
    mu_new<-(sig1^2*mu0*(1-lambda) + sig0^2*mu1*lambda) / (sig1^2*(1-lambda) + sig0^2*lambda)
    return(mu_new)
}

# sigma update rule for normals:
sig_update<-function(lambda, T) {
    sig_new<-sqrt( T*(sig0^2*sig1^2) / (sig1^2*(1-lambda) + sig0^2*lambda) )
    return(sig_new)
}

# make a function that will take lambda, T and a draw to generate an unnormalized density
unnorm_dens<-function(x, lambda, T) {
    mu_new<-mu_update(lambda, T)
    sigma_new<-sig_update(lambda, T)
    q_dens<-exp(-(x-mu_new)^2/(2*sigma_new^2))
    return(q_dens)
}

##############
# movement options
##############

# define degrees of freedom and maximum displacement for each

dof<-c('TRUE', temp_dof)
names(dof)<-c('lambda', 'temp')

dof_ranges<-list()
dof_ranges[1]<-list(0:move_jump)
if (dof[2]) { dof_ranges[2]<-list(-move_jump : move_jump)
    }  else { dof_ranges[2]<-list(NA) }
names(dof_ranges)<-names(dof)

nacheck<-function(x) { return(any(is.na(x))) }

get_all_moves<-function() {
    grid_full<-expand.grid(dof_ranges)
    grid_strip<-grid_full[ , !apply(grid_full, 2, nacheck), drop=F] # remove NA columns that mess with absolute value
    grid_trim<-grid_strip[rowSums(abs(grid_strip))<=move_jump, , drop=F] # filter out actions with too many moves
    grid_trim2<-grid_trim[rowSums(abs(grid_trim))>0, , drop=F] # filter out null moves
    ind<-as.numeric(rownames(grid_trim2)) # pull out good rows
    grid_out<-grid_full[ind, ]
    return(grid_out)
}

move_universe<-get_all_moves()

equivcheck<-function(list, point) { isTRUE(all.equal(list, point)) }

# movement primitive
move_result<-function(move, state) {

    # compute new state based on original state and set of changes from move
    new_pos<-rep(NA, length.out=length(state))

    for (i in 1:length(move)) { # for each dimension
        if (!is.na(move[i])) { # if that dimension is varying
            dimension<-point_lists[[i]] # get possible values
            curr_ind<-which(sapply(dimension, equivcheck, state[i]))
            new_ind<-curr_ind + as.numeric(move[i])
            if (new_ind<=0) { new_pos[i]<-NA # if out of bounds, call it NA
            } else if (new_ind>length(dimension)) { new_pos[i]<-NA 
            } else { new_pos[i]<-dimension[new_ind] }
        } else { new_pos[i]<-state[i] } # if move is NA, then retain position
    } # note: output may have NAs in position - this is fine. NA moves are pruned by valid move checker
    return(new_pos)
}

is_valid_move<-function(move, state) { 
    pos<-paste(move_result(move, state), collapse='_')
    if (pos %in% indexer) { return(TRUE) 
    } else { return(FALSE) }
}

# create adjacency lists
adj_list<-vector('list', length(indexer))
names(adj_list)<-indexer
for (i in 1:length(point_list)) {
    move_bool<-apply(move_universe, 1, is_valid_move, point_list[[i]])
    valid_moves<-move_universe[move_bool,]
    neighbors<-t(apply(valid_moves, 1, move_result, point_list[[i]]))

    # if first entry of point list is lambdamax, then only allow neighbors that are above/below such that it goes towards the target
    if (point_list[[i]][1]==lambda_max) {
        rm_ind<-vector()
        if (point_list[[i]][2]>state_target[2]) {
            rm_ind<-which(neighbors[,2]>point_list[[i]][2])
        } else if (point_list[[i]][2]<state_target[2]) {
            rm_ind<-which(neighbors[,2]<point_list[[i]][2])
        }
        if (length(rm_ind)!=0) { neighbors<-neighbors[-rm_ind, , drop=FALSE] }
    }
    adj_list[[i]]<-apply(neighbors, 1, paste, collapse='_')
}


##############
# map setups
##############

q_map<-vector('list', length(indexer))
names(q_map)<-indexer
for (i in 1:length(point_list)) {

    # for each point in the graph, create an entry to the Q map
    # each entry is a matrix denoting the legal actions/next states and the Q value associated with that action
    q_mat<-array(NA, dim=c(length(adj_list[[i]]), 2))
    q_mat[,1]<-adj_list[[i]]
    q_mat[,2]<-0 # arbitrary large value to bias path towards explored edges towards end of algorithm (decision reversed, see below)
    q_map[[i]]<-q_mat
}

# consider initializing q to 0. 
# need special handling for maps at target states - no moves are allowed, and the q score for adjacent states should just be the discounted reward

r_map<-vector('list', length(indexer))
names(r_map)<-indexer
for (i in 1:length(point_list)) {

    # for each point in the graph create an entry to the Rewards map
    # each entry is a list indexed by the next states, where each list entry is a matrix showing guesses for the reward and their variances (if available)
    r_list<-vector('list', length(adj_list[[i]]))
    names(r_list)<-adj_list[[i]]
    for (j in 1:length(adj_list[[i]])) {
        r_list[[j]]<-array(NA, dim=c(1,2)) # initialize matrix with NAs for correct dimension. use na.rm option to remove them from future calculations
    }
    r_map[[i]]<-r_list
}

get_valid_moves<-function(state, path) {

    # get all adjacent states
    curr_ind<-which(indexer==state)
    full_neighbors<-adj_list[[curr_ind]]

    # remove states that are present in the path to date
    trim_neighbors<-setdiff(full_neighbors, path)

    return(trim_neighbors)
}

#############
# sampling helpers
#############

sis_transition<-function(draw, distr_indexer) {
    A<-draw
    for (i in 1:nummetrostep) {
        for (j in 1:length(metro_spread)) {
            B<-metro_func(A, distr_indexer, metro_spread[j])
            A<-B
        }
    }
    return(B)
}

metro_func<-function(draw, distr_indexer, spread) {
    trial_draw<-rnorm(1, mean=draw, sd=spread)
    distr_state<-as.numeric(unlist(strsplit(distr_indexer, split='_')))
    curr_dens<-unnorm_dens(draw, distr_state[1], distr_state[2])
    trial_dens<-unnorm_dens(trial_draw, distr_state[1], distr_state[2])
    accept_prob<-min(1, trial_dens/curr_dens)
    prob_draw<-runif(n=1, min=0, max=1)
    if (prob_draw<=accept_prob) { return(trial_draw) 
    } else { return(draw) }
}

get_incr_weight<-function(draw, indexer_prev, indexer_curr) {
    distr_prev<-as.numeric(unlist(strsplit(indexer_prev, split='_')))
    distr_curr<-as.numeric(unlist(strsplit(indexer_curr, split='_')))
    curr_dens<-unnorm_dens(draw, distr_curr[1], distr_curr[2])
    prev_dens<-unnorm_dens(draw, distr_prev[1], distr_prev[2])
    incr_weights<-curr_dens/prev_dens
    return(incr_weights)
}

calc_norm_ratio<-function(norm_weights, incr_weights) {
    return(sum(norm_weights*incr_weights))
}

boot_var_ratio<-function(norm_weights, incr_weights) {

    # sample with replacement to get new set
    boot_index<-sample(1:length(norm_weights), size=length(norm_weights), replace=T)
    boot_norm<-norm_weights[boot_index]/sum(norm_weights[boot_index])
    boot_incr<-incr_weights[boot_index]
    boot_ratio<-calc_norm_ratio(boot_norm, boot_incr)
}