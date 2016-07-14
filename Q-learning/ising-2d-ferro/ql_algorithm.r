# ql_algorithm.r

# to be called by ql_master.r
# for use with the 2D Ising ferromagnet example
# define Q-learning algorithm and helpers. adapted to minimize cost instead of maximize reward.

ql_full_history<-function(q_map, r_map, alpha, gamma, epsilon_init, epsilon_tau) {

    maps<-list(q_map, r_map)
    delta_epsilon<-1/epsilon_tau
    path_solns<-vector('list', length=num_episode)
    path_ratios<-rep(NA, length.out=num_episode)
    path_vars<-rep(NA, length.out=num_episode)

    for (i in 1:epsilon_tau) {

        # in first phase, epsilon varies predictably from fully random to fully greedy

        # print progress
        print(i)

        # run one episode and update maps
        maps<-ql_episode(maps[[1]], maps[[2]], epsilon_init+(i-1)*delta_epsilon)

        # get summary output - solution paths, ratios, variances...
        path_soln_temp<-ql_path_soln(maps[[1]])
        path_ratio_temp<-ql_path_ratio(maps[[2]], path_soln_temp)
        path_var_temp<-ql_path_var(maps[[2]], path_soln_temp)
        path_soln_mat<-t(sapply(path_soln_temp, split_func))
        if (!any(is.na(path_soln_mat))) { colnames(path_soln_mat)<-c('beta', 'magnetic_field') }
        path_solns[[i]]<-path_soln_mat
        path_ratios[i]<-path_ratio_temp
        path_vars[i]<-path_var_temp
    }

    for (i in (epsilon_tau+1):num_episode) {

        # in the second phase, epsilon is fixed on the fully greedy value. exploration is over.

        # print progress
        print(i)

        # run episode
        maps<-ql_episode(maps[[1]], maps[[2]], 1)

        # get summary output
        path_soln_temp<-ql_path_soln(maps[[1]])
        path_ratio_temp<-ql_path_ratio(maps[[2]], path_soln_temp)
        path_var_temp<-ql_path_var(maps[[2]], path_soln_temp)
        path_soln_mat<-t(sapply(path_soln_temp, split_func))
        if (!any(is.na(path_soln_mat))) { colnames(path_soln_mat)<-c('beta', 'magnetic_field') }
        path_solns[[i]]<-path_soln_mat
        path_ratios[i]<-path_ratio_temp
        path_vars[i]<-path_var_temp
    }

    # return full history of summary output
    return(list(path_solns, path_ratios, path_vars))
}

ql_episode<-function(q_map, r_map, epsilon) {

    # ql_episode subroutine represents a single episode of learning
    # each wave explores a single path using one bundle of samples.

    # initialize this episode's search agent
    # agent holds info on current path and draws
    # at end, reprocess agent draws to determine ratios and vars, then update maps

    ##########
    # step 1: create path and collect forward trajectory draws
    ##########

    agent_path<-vector()
    agent_path[1]<-indexer_init
    iter_dummy<-0

    # create storage for forward draws and weights
    draws_F<-array(NA, dim=c(num_traj, lattice_size^2, 1))
    weights_F<-array(NA, dim=c(num_traj, 1))

    # initialize values
    draws_F[ , , 1]<-generate_ising_hightemp(num_traj, lattice_size)
    weights_F[,1]<-1/num_traj

    while(!identical(agent_path[length(agent_path)], indexer_target)) { # until the agent has reached the target, continue the search

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

        # propagate particles
        new_draws<-t(apply(draws_F[ , , iter_dummy], 1, ising_transition, indexer_next))
        draws_F<-abind(draws_F, new_draws, along=3)

        # calculate particle weights
        p_weight<-get_p_weight(draws_F[ , , iter_dummy], weights_F[,iter_dummy], agent_path[iter_dummy], agent_path[iter_dummy+1])
        weights_F<-cbind(weights_F, p_weight)

        # continue particle propagation until path is complete
    }

    # propagate particles backwards
    draws_R<-array(NA, dim=dim(draws_F))
    weights_R<-array(NA, dim=dim(weights_F))

    ######
    # initialize randomly to mostly spin up or down
    ######
    draws_R[ , , dim(draws_F)[3]]<-generate_ising_lowtemp(num_traj, lattice_size, sample(c(-1, 1), size=1))
    weights_R[, ncol(weights_R)]<-1/num_traj

    for (i in (ncol(weights_R)-1):1) {
        draws_R[ , , i]<-t(apply(draws_R[, , i+1], 1, ising_transition, agent_path[i]))
        weights_R[,i]<-get_p_weight(draws_R[ , , i+1], weights_R[, i+1], agent_path[i+1], agent_path[i])
    }

    # process forward and reverse weighted trajectories to obtain ratio and var(ratio) estimates
    #####
    # NOTE: IF THIS IS FOR ONLY THE POSITIVE MAGNETIZATION CASE, WE ARE LIKELY MISSING OUT ON 1/2 OF THE TARGET PARTITION FUNCTION
    #####
    r_map_updates<-process_trajectories(draws_F, draws_R, weights_F, weights_R, agent_path)

    # update maps
    for (i in 1:ncol(r_map_updates)) {

        # get current state and next state indexers and indices
        state_current<-agent_path[i]
        state_next<-agent_path[i+1]
        index_current<-which(indexer==state_current)
        index_next<-which(names(r_map[[index_current]])==state_next)

        # add r map entries
        r_map[[index_current]][[index_next]]<-rbind(r_map[[index_current]][[index_next]], c(r_map_updates[,i], 1))

        # update q map
        q_entry<-as.numeric(q_map[[index_current]][index_next, 2])
        r_entry<-mean(r_map[[index_current]][[index_next]][,2] , na.rm=TRUE)
        q_prime_ind<-which(indexer==state_next)
        q_prime_val<-min(as.numeric(q_map[[q_prime_ind]][,2]))
        q_map[[index_current]][index_next, 2]<-q_score(q_entry, r_entry, q_prime_val)
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
    if (q==Inf) { # if the previous q value was undefined (Inf)
        if (q_prime==Inf) { # if the next step q is also undefined (Inf)
            new_qsc<-r+r # set alpha=1 for instant learning, and use current r as best guess for next q
        } else {
            new_qsc<-r+gamma*q_prime # otherwise just set alpha=1 for instant learning
        }
    } else if (q_prime==Inf) { # if only next step q is undefined
        new_qsc<- q + alpha*(r+r-q)
    } else { # if neither is undefined then just go on as normal
        new_qsc<-q + alpha*(r + gamma*q_prime - q)
    }
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

ql_path_ratio<-function(r_map, path_soln) {
    ratio<-1
    if (any(is.na(path_soln))) { ratio<-NA
    } else {
        for (i in 1:(length(path_soln)-1)) {
            prev_ind<-which(indexer==path_soln[[i]])
            ratio_ind<-which(names(r_map[[prev_ind]])==path_soln[i+1])
            step_ratio<-weighted.mean(r_map[[prev_ind]][[ratio_ind]][,1], 1/r_map[[prev_ind]][[ratio_ind]][,2], na.rm=TRUE)
            ratio<-ratio*step_ratio
        }
    }
    return(ratio)
}

ql_path_var<-function(r_map, path_soln) {
    cost<-0
    if (any(is.na(path_soln))) { cost<-NA
    } else {
        for (i in 1:(length(path_soln)-1)) {
            prev_ind<-which(indexer==path_soln[[i]])
            cost_ind<-which(names(r_map[[prev_ind]])==path_soln[i+1])
            step_cost<-weighted.mean(r_map[[prev_ind]][[cost_ind]][,2], r_map[[prev_ind]][[cost_ind]][,3], na.rm=TRUE)
            cost<-cost+step_cost
        }
    }
    return(cost)
}

ql_path_var_scaled<-function(r_map, path_soln) {
    cost<-0
    if (any(is.na(path_soln))) { cost<-NA
    } else {
        for (i in 1:(length(path_soln)-1)) {
            prev_ind<-which(indexer==path_soln[[i]])
            cost_ind<-which(names(r_map[[prev_ind]])==path_soln[i+1])
            step_cost<-weighted.mean(r_map[[prev_ind]][[cost_ind]][,2], r_map[[prev_ind]][[cost_ind]][,3], na.rm=TRUE)
            step_traversals<-nrow(r_map[[prev_ind]][[cost_ind]])-1 # remove one for the NA row used for initialization
            cost<-cost+step_cost/step_traversals
        }
    }
    return(cost)
}

split_func<-function(x) { 
    if(any(is.na(x))) { out<-NA
    } else { out<-as.numeric(unlist(strsplit(x, split='_'))) }
    return(out)
}




