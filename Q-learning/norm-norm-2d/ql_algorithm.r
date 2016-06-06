# ql_algorithm.r

# to be called by ql_master.r
# define Q-learning algorithm. adapted to minimize cost instead of maximize reward

# add a piece to ql_search that checks for convergence
# crtieria: steady opt path
# criteria: steady path cost
# note: need function to return opt path (handle incomplete ones)
# note: need to store prev solns for comparison
# note: need new handling for epsilon - no longer based on a static schedule since we don't know total num of waves

ql_search<-function(q_map, r_map, alpha, gamma, epsilon_init, min_episode, conv_tol) {

    # run QL episodes until path and path cost converge
    maps<-list(q_map, r_map)
    
    # run the minimum number of searches
    for (ql_iter in 1:min_episode) {
        print(ql_iter)
        maps<-ql_episode(maps[[1]], maps[[2]], (i-1)*epsilon_delta)
    }

    # calculate the current optimal path and cost
    prev_path<-ql_path_soln(maps[[1]])
    prev_cost<-ql_path_cost(maps[[2]], prev_path)    

    # set convergence flags to false
    path_conv<-FALSE
    cost_conv<-FALSE
    epsilon_curr<-0.5

    while (!(path_conv & cost_conv)) { # until both are converged, keep going
        ql_iter<-ql_iter+1
        print(ql_iter)

        # run an additional round and update info
        maps<-ql_episode(maps[[1]], maps[[2]], epsilon_curr)
        curr_path<-ql_path_soln(maps[[1]])
        curr_cost<-ql_path_cost(maps[[2]], curr_path)

        # check for convergence
        path_conv<-identical(prev_path, curr_path)
        cost_delta<-abs((curr_cost-prev_cost)/prev_cost)
        cost_conv<-(cost_delta<conv_tol)

        # update info for next round for searching
        prev_path<-curr_path
        prev_cost<-curr_cost
        epsilon_curr<-epsilon_update(epsilon_curr, path_conv, cost_delta)
    }
    return(maps)
}

## ql_episode subroutine represents a single episode of learning
# each wave explores a single path using one bundle of samples
# inputs: 
# r_map : list by indexer (one entry per node), within list, place list by adjacent nodes, place n x 2 matrix for varBAR and varvarBAR for each bundle's evaluation
# q_map : list by indexer (one entry per node), each list entry  is vector of q values to neighbors

ql_episode<-function(q_map, r_map, epsilon) {

    # initialize this episode's search agent
    # agent holds info on current path and draws
    # at end, reprocess agent draws to determine var of each edge
    # report vars as rewards and update q map
    agent_path<-vector()

    agent_path[1]<-indexer_init
    iter_dummy<-0

    while (!identical(agent_path[length(agent_path)], indexer_target)) { # until the agent has reached the target, continue the search

        iter_dummy<-iter_dummy+1

        # set current state
        state_current<-agent_path[length(agent_path)]
        curr_index<-which(indexer==state_current)

        # get valid moves and q scores for those moves
        full_q<-q_map[[curr_index]]
        trim_index<-which(full_q[,1] %in% agent_path)
        if (length(trim_index)!=0) { trim_q<-full_q[-trim_index, , drop=FALSE] 
        } else { trim_q<-full_q }

        # select between those moves based on q score policy
        indexer_next<-q_choice(trim_q, epsilon)
        next_index<-which(names(r_map[[curr_index]])==indexer_next)

        # add choice to path
        agent_path[iter_dummy+1]<-indexer_next

        # draw new edge reward from distribution
        reward_entry<-c(rnorm(1, r_distns[[curr_index]][next_index,1], r_distns[[curr_index]][next_index,2]),1)

        # update r_map
        r_map[[curr_index]][[next_index]]<-rbind(r_map[[curr_index]][[next_index]], reward_entry)

        # update q_map
        q_entry<-as.numeric(q_map[[curr_index]][next_index, 2])
        r_entry<-mean(r_map[[curr_index]][[next_index]][,1] , na.rm=TRUE)
        q_prime_ind<-which(indexer==indexer_next)
        q_prime_val<-min(as.numeric(q_map[[q_prime_ind]][,2]))
        q_map[[curr_index]][next_index, 2]<-q_score(q_entry, r_entry, q_prime_val)
    }
    return(list(q_map, r_map))
}

# epsilon-decision making policy
# with probability epsilon, select next state based on best q
# with prob 1-epsilon, select randomly
q_choice<-function(map, e) {
    rand_val<-runif(1)
    if (rand_val <= e) { # choose best q
        return(map[which.min(map[,2]),1])
    } else { return(map[sample(1:nrow(map),1),1]) }
}

q_score<-function(q, r, q_prime) {
    new_qsc<-q + alpha*(r + gamma*q_prime - q)
    return(new_qsc)
}

ql_path_soln<-function(q_map) {
    path_soln<-list(indexer_init)

    # until final state in path is target state, keep searching unless have break
    while (path_soln[[length(path_soln)]] != indexer_target) {
        prev_ind<-which(indexer==path_soln[[length(path_soln)]])
        prev_q<-q_map[[prev_ind]]
        next_ind<-which.min(prev_q[,2])
        next_path<-prev_q[next_ind,1]

        if (next_path %in% path_soln) { # if path contains loop, break and return no path
            path_soln<-NA
            break
        } else {
            path_soln[[length(path_soln)+1]]<-next_path
        }
    }
    return(path_soln)
}

# potentially replace with a score-convergence based on the Q map?
ql_path_cost<-function(r_map, path_soln) {
    cost<-0
    if (any(is.na(path_soln))) { cost<-NA 
    } else {
        for (i in 1:(length(path_soln)-1)) {
            prev_ind<-which(indexer==path_soln[[i]])
            cost_ind<-which(names(r_map[[prev_ind]])==path_soln[i+1])
            step_cost<-mean(r_map[[prev_ind]][[cost_ind]][,1] , na.rm=TRUE)
            cost<-cost+step_cost
        }
    }   
    return(cost)
}

epsilon_update<-function(epsilon_curr, path_conv, cost_delta) {
    if (path_conv) { # if path is converged
        if (cost_delta>0.25) { # but value is really noisy, don't change randomness
            epsilon_new<-epsilon_curr
        } else { epsilon_new<-epsilon_curr+0.025 } # if value is reasonably stable, reduce randomness
    } else { #if path isn't converged
        if (cost_delta<0.1) { # but the values between competing paths are reasonably close, reduce randomness
            epsilon_new<-epsilon_curr+0.025
        } else { epsilon_new<-epsilon_curr } 
    }
    return(epsilon_new)
}

#### for testing only
ql_path_cost2<-function(r_map, path_soln) {
    cost<-list()
    if (any(is.na(path_soln))) { cost<-NA 
    } else {
        for (i in 1:(length(path_soln)-1)) {
            prev_ind<-which(indexer==path_soln[[i]])
            cost_ind<-which(names(r_map[[prev_ind]])==path_soln[i+1])
            step_cost<-mean(r_map[[prev_ind]][[cost_ind]][,1] , na.rm=TRUE)
            cost[[i]]<-step_cost
        }
    }   
    return(cost)
}