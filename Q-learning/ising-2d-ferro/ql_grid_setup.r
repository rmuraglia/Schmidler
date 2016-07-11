# ql_grid_setup.r

# set up grid for Q-learning for Ising model in beta-magn parameter space.
# to be called by ql_master.r - has depends on params defined in this file

##################
# set up grid
##################

beta_points<-seq(from=beta_min, to=beta_max, length.out=beta_numpoints)
# magn_points<-seq(from=-magn_max_abs, to=magn_max_abs, length.out=magn_numpoints) # external magnetic field intensities are symmetric about zero.
magn_points<-seq(from=magn_min, to=magn_max_abs, length.out=magn_numpoints)

# list of length 2 containing state space param values
point_lists<-list(beta_points, magn_points)
names(point_lists)<-c('beta', 'magnetic_field')

# matrix containing combinations of state space param values
point_grid<-expand.grid(point_lists)

# list of combinations of param values
point_list<-as.list(as.data.frame(t(point_grid)))

# unique string identifiers for each param combination
indexer<-apply(point_grid, 1, paste, collapse='_')

# set target magnetization and key points
state_init<-c(beta_min, magn_min)
state_target<-c(beta_max, magn_min)
indexer_init<-paste(state_init, collapse='_')
indexer_target<-paste(state_target, collapse='_')

##############
# movement options
##############

# define degrees of freedom and maximum displacement for each

dof<-c(TRUE, magn_dof)
names(dof)<-c('beta', 'magnetic_field')

dof_ranges<-list()
dof_ranges[1]<-list(0:move_jump) # for beta, disallow left moves. must progress toward target beta
if (dof[2]) { dof_ranges[2]<-list(-move_jump : move_jump) # if magn is allowed to vary, it may go up or down.
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

    # if first entry of point list is betamax, then only allow neighbors that are above/below such that it goes towards the target
    if (magn_dof && point_list[[i]][1]==beta_max) {
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
    # q_mat[,2]<-1 # arbitrary large value to bias path towards explored edges towards end of algorithm (decision reversed, see below)
    q_mat[,2]<-Inf # try infinite intialized value for when you have no idea of the scales
    q_map[[i]]<-q_mat
}

q_map[[which(names(q_map)==indexer_target)]][,2]<-0

r_map<-vector('list', length(indexer))
names(r_map)<-indexer
for (i in 1:length(point_list)) {

    # for each point in the graph create an entry to the Rewards map
    # each entry is a list indexed by the next states, where each list entry is a matrix showing guesses for the ratio, var(ratio) and var(var[ratio])  (if available)
    r_list<-vector('list', length(adj_list[[i]]))
    names(r_list)<-adj_list[[i]]
    for (j in 1:length(adj_list[[i]])) {
        r_list[[j]]<-array(NA, dim=c(1,3)) # initialize matrix with NAs for correct dimension. use na.rm option to remove them from future calculations
    }
    r_map[[i]]<-r_list
}

