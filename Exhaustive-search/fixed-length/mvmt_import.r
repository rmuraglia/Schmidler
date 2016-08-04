# mvmt_import.r

# degrees of freedom: lambda can always change, but sometimes restrict temperature
dof<-c('TRUE', temp_dof) 
names(dof)<-c('lambda', 'temperature')

# get maximum range of movement options in each dimension
dof_ranges<-list()
for (i in dof) {
    if (i) { dof_ranges<-append(dof_ranges, list(-move_jump: move_jump))
    } else { dof_ranges<-append(dof_ranges, list(NA)) }
}
names(dof_ranges)<-names(dof)

# helper function to trim locked dimensions
nacheck<-function(x) { return(any(is.na(x))) }

# get all possible moves given jump coefficient and degrees of freedom, ignoring position on grid
get_all_moves<-function() {
    grid_full<-expand.grid(dof_ranges) # complete listing of moves
    grid_strip<-grid_full[ , !apply(grid_full, 2, nacheck), drop=F] # remove locked dimension columns
    grid_trim1<-grid_strip[rowSums(abs(grid_strip))<=move_jump, , drop=F] # filter out actions with too many moves
    grid_trim2<-grid_trim1[rowSums(abs(grid_trim1))>0, , drop=F] # filter out null moves
    ind<-as.numeric(rownames(grid_trim2)) # get indices of surviving moves
    grid_out<-grid_full[ind,]
    return(grid_out)
}

move_universe<-get_all_moves()

# helper function to check equivalency
equivcheck<-function(list, point) { isTRUE(all.equal(list, point)) }

# movement primitive
move_result<-function(action, curr_state) {
    # compute new state based on curr_state and action
    new_state<-rep(NA, length.out=length(curr_state))

    for (i in 1:length(action)) { # for each dimension
        if (!is.na(action[i])) { # if dimension is not locked
            dimension<-point_lists[[i]] # get potential values for dimension
            curr_ind<-which(sapply(dimension, equivcheck, curr_state[i])) # get current state's index for that dimension
            new_ind<-curr_ind + as.numeric(action[i]) # get new index based on action
            if (new_ind<=0) { new_state[i]<-NA # if out of bounds, call it NA
            } else if (new_ind>length(dimension)) { new_state[i]<-NA 
            } else { new_state[i]<-dimension[new_ind] } # otherwise get new position
        } else { new_state[i]<-curr_state[i] } # if dimension IS locked, then stay where you are
    } # note: output may have NAs in position. will get pruned later
    return(new_state)
}

# check if a given move is valid
is_valid_move<-function(action, curr_state) {
    pos<-paste(move_result(action, curr_state), collapse='_')
    if (pos %in% indexer) { return(TRUE)
    } else { return(FALSE) }
}