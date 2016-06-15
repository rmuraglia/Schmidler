# ql_grid_setup.r
# maze grid version - doesn't need all the SIS and distribution portions

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
    q_mat[,2]<-100 # arbitrary large value to bias path towards explored edges towards end of algorithm (decision reversed, see below)
    q_map[[i]]<-q_mat
}

q_map[[which(names(q_map)==indexer_target)]][,2]<-0

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

###################
# user defined reward map
###################

# manually create a reward map lookup, so we can construct an example with known paths.
# each reward is normally distributed with varying centers and scales
# each traversal of an edge will draw a random reward from the corresponding distribution.

r_distns<-vector('list', length(indexer))
names(r_distns)<-indexer

# for each state, create an entry to the list, detailing cost and variance (noise) for each exiting edge
# 0, 0.5
k<-1
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(10,1)
r_d_temp[2,]<-c(1,1)
r_distns[[k]]<-r_d_temp

# 0.2, 0.5
k<-2
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(50,1)
r_d_temp[2,]<-c(2.5,1)
r_distns[[k]]<-r_d_temp

# 0.4, 0.5
k<-3
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(10,1)
r_d_temp[2,]<-c(1,1)
r_distns[[k]]<-r_d_temp

# 0.6, 0.5
k<-4
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(50,1)
r_d_temp[2,]<-c(1,1)
r_distns[[k]]<-r_d_temp

# 0.8, 0.5
k<-5
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(10,1)
r_d_temp[2,]<-c(2.5,1)
r_distns[[k]]<-r_d_temp

# 1, 0.5
k<-6
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(1,1)
r_distns[[k]]<-r_d_temp

# 0, 1.375
k<-7
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(1,1)
r_d_temp[2,]<-c(10,1)
r_d_temp[3,]<-c(1,1)
r_distns[[k]]<-r_d_temp

# 0.2, 1.375
k<-8
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(2.5,1)
r_d_temp[2,]<-c(50,1)
r_d_temp[3,]<-c(2.5,1)
r_distns[[k]]<-r_d_temp

# 0.4, 1.375
k<-9
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(1,1)
r_d_temp[2,]<-c(1,1)
r_d_temp[3,]<-c(1,1)
r_distns[[k]]<-r_d_temp

# 0.6, 1.375
k<-10
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(1,1)
r_d_temp[2,]<-c(50,1)
r_d_temp[3,]<-c(1,1)
r_distns[[k]]<-r_d_temp

# 0.8, 1.375
k<-11
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(2.5,1)
r_d_temp[2,]<-c(10,1)
r_d_temp[3,]<-c(2.5,1)
r_distns[[k]]<-r_d_temp

# 1, 1.375
k<-12
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(1,1)
r_d_temp[2,]<-c(1,1)
r_distns[[k]]<-r_d_temp

# 0, 2.25
k<-13
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(1,1)
r_d_temp[2,]<-c(10,1)
r_d_temp[3,]<-c(1,1)
r_distns[[k]]<-r_d_temp

# 0.2, 2.25
k<-14
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(2.5,1)
r_d_temp[2,]<-c(1,1)
r_d_temp[3,]<-c(2.5,1)
r_distns[[k]]<-r_d_temp

# 0.4, 2.25
k<-15
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(1,1)
r_d_temp[2,]<-c(10,1)
r_d_temp[3,]<-c(1,1)
r_distns[[k]]<-r_d_temp

# 0.6, 2.25
k<-16
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(1,1)
r_d_temp[2,]<-c(50,1)
r_d_temp[3,]<-c(1,1)
r_distns[[k]]<-r_d_temp

# 0.8, 2.25
k<-17
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(2.5,1)
r_d_temp[2,]<-c(10,1)
r_d_temp[3,]<-c(2.5,1)
r_distns[[k]]<-r_d_temp

# 1, 2.25
k<-18
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(1,1)
r_distns[[k]]<-r_d_temp

# 0, 3.125
k<-19
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(1,1)
r_d_temp[2,]<-c(10,1)
r_d_temp[3,]<-c(1,1)
r_distns[[k]]<-r_d_temp

# 0.2, 3.125
k<-20
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(2.5,1)
r_d_temp[2,]<-c(50,1)
r_d_temp[3,]<-c(2.5,1)
r_distns[[k]]<-r_d_temp

# 0.4, 3.125
k<-21
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(1,1)
r_d_temp[2,]<-c(10,1)
r_d_temp[3,]<-c(1,1)
r_distns[[k]]<-r_d_temp

# 0.6, 3.125
k<-22
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(1,1)
r_d_temp[2,]<-c(1,1)
r_d_temp[3,]<-c(1,1)
r_distns[[k]]<-r_d_temp

# 0.8, 3.125
k<-23
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(2.5,1)
r_d_temp[2,]<-c(10,1)
r_d_temp[3,]<-c(2.5,1)
r_distns[[k]]<-r_d_temp

# 1, 3.125
k<-24
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(1,1)
r_distns[[k]]<-r_d_temp

# 0, 4
k<-25
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(1,1)
r_d_temp[2,]<-c(10,1)
r_distns[[k]]<-r_d_temp

# 0.2, 4
k<-26
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(2.5,1)
r_d_temp[2,]<-c(50,1)
r_distns[[k]]<-r_d_temp

# 0.4, 4
k<-27
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(1,1)
r_d_temp[2,]<-c(10,1)
r_distns[[k]]<-r_d_temp

# 0.6, 4
k<-28
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(1,1)
r_d_temp[2,]<-c(50,1)
r_distns[[k]]<-r_d_temp

# 0.8, 4
k<-29
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(2.5,1)
r_d_temp[2,]<-c(10,1)
r_distns[[k]]<-r_d_temp

# 1, 4
k<-30
r_d_temp<-array(NA, dim=dim(q_map[[k]]))
r_d_temp[1,]<-c(1,1)
r_distns[[k]]<-r_d_temp