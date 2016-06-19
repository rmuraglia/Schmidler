# ql_algorithm.r

# to be called by ql_master.r
# define Q-learning algorithm. adapted to minimize cost instead of maximize reward

ql_full_history<-function(q_map, r_map, alpha, gamma, epsilon_init, epsilon_tau) {

    maps<-list(q_map, r_map)
    delta_epsilon<-1/epsilon_tau
    path_solns<-vector('list', length=num_episode)
    path_ratios<-rep(NA, length.out=num_episode)
    path_vars<-rep(NA, length.out=num_episode)

    for (i in 1:epsilon_tau) {
        print(i)
        maps<-ql_episode(maps[[1]], maps[[2]], epsilon_init+(i-1)*delta_epsilon)
        path_soln_temp<-ql_path_soln(maps[[1]])
        path_ratio_temp<-ql_path_ratio(maps[[2]], path_soln_temp)
        path_var_temp<-ql_path_var(maps[[2]], path_soln_temp)
        path_soln_mat<-t(sapply(path_soln_temp, split_func))
        if (!any(is.na(path_soln_mat))) { colnames(path_soln_mat)<-c('lambda', 'temperature') }
        path_solns[[i]]<-path_soln_mat
        path_ratios[i]<-path_ratio_temp
        path_vars[i]<-path_var_temp
    }
    for (i in (epsilon_tau+1):num_episode) {
        print(i)
        maps<-ql_episode(maps[[1]], maps[[2]], 1)
        path_soln_temp<-ql_path_soln(maps[[1]])
        path_ratio_temp<-ql_path_ratio(maps[[2]], path_soln_temp)
        path_var_temp<-ql_path_var(maps[[2]], path_soln_temp)
        path_soln_mat<-t(sapply(path_soln_temp, split_func))
        if (!any(is.na(path_soln_mat))) { colnames(path_soln_mat)<-c('lambda', 'temperature') }
        path_solns[[i]]<-path_soln_mat
        path_ratios[i]<-path_ratio_temp
        path_vars[i]<-path_var_temp
    }
    return(list(path_solns, path_ratios, path_vars))
}

ql_three_chain<-function(q_map, r_map, alpha, gamma, min_episode, conv_tol) {

    # initialize three map sets
    maps1<-list(q_map, r_map)
    maps2<-maps1
    maps3<-maps1
    epsilons<-c(0,0,0)

    # run minimum number of random searches
    for (ql_iter in 1:min_episode) {
        print(ql_iter)
        maps1<-ql_episode(maps1[[1]], maps1[[2]], 0)
        maps2<-ql_episode(maps2[[1]], maps2[[2]], 0)
        maps3<-ql_episode(maps3[[1]], maps3[[2]], 0)
    }

    # get solution paths, ratio estimates and variance estimates
    soln1<-ql_path_soln(maps1[[1]])
    soln2<-ql_path_soln(maps2[[1]])
    soln3<-ql_path_soln(maps3[[1]])
    ratio1<-ql_path_ratio(maps1[[2]], soln1)
    ratio2<-ql_path_ratio(maps2[[2]], soln2)
    ratio3<-ql_path_ratio(maps3[[2]], soln3)
    var1<-ql_path_var(maps1[[2]], soln1)
    var2<-ql_path_var(maps2[[2]], soln2)
    var3<-ql_path_var(maps3[[2]], soln3)

    ## if desired, add convergence criteria based on path identity
    # path_match<-c(identical(soln1, soln2), identical(soln1, soln3), identical(soln2, soln3))
    # path_conv<-all(path_match)
    path_conv<-TRUE

    ## convergence criteria based on ratio estimate similarity
    # get standard deviations and ratio intervals
    sd1<-sqrt(var1*ratio1^2)
    sd2<-sqrt(var2*ratio2^2)
    sd3<-sqrt(var3*ratio3^2)
    interval1<-c(ratio1-sd1*conv_tol, ratio1+sd1*conv_tol)
    interval2<-c(ratio2-sd2*conv_tol, ratio2+sd2*conv_tol)
    interval3<-c(ratio3-sd3*conv_tol, ratio3+sd3*conv_tol)

    # determine if intervals are overlapping
    overlap1<-overlap_check(interval1[1], interval1[2], interval2[1], interval2[2])
    overlap2<-overlap_check(interval1[1], interval1[2], interval3[1], interval3[2])
    overlap3<-overlap_check(interval2[1], interval2[2], interval3[1], interval3[2])
    cost_match<-c(overlap1, overlap2, overlap3)
    cost_conv<-all(cost_match)

    # storage for cost_matches
    cost_archive<-array(NA, dim=c((max_episode-min_episode), length(cost_match)))

    # set epsilons to a mostly random value
    epsilons<-c(0.3, 0.3, 0.3)

    while(!(path_conv & cost_conv)) {
        ql_iter<-ql_iter+1
        print(ql_iter)

        # update epsilons based on cost match only
        epsilons<-epsilon_update_cost(epsilons, cost_match)

        # run next search episode
        maps1<-ql_episode(maps1[[1]], maps1[[2]], epsilons[1])
        maps2<-ql_episode(maps2[[1]], maps2[[2]], epsilons[2])
        maps3<-ql_episode(maps3[[1]], maps3[[2]], epsilons[3])

        # get new solution paths, ratio estimates and variance estimates
        soln1<-ql_path_soln(maps1[[1]])
        soln2<-ql_path_soln(maps2[[1]])
        soln3<-ql_path_soln(maps3[[1]])
        ratio1<-ql_path_ratio(maps1[[2]], soln1)
        ratio2<-ql_path_ratio(maps2[[2]], soln2)
        ratio3<-ql_path_ratio(maps3[[2]], soln3)
        var1<-ql_path_var(maps1[[2]], soln1)
        var2<-ql_path_var(maps2[[2]], soln2)
        var3<-ql_path_var(maps3[[2]], soln3)

        # update convergence criteria
        # path_match<-c(identical(soln1, soln2), identical(soln1, soln3), identical(soln2, soln3))
        # path_conv<-all(path_match)
        sd1<-sqrt(var1*ratio1^2)
        sd2<-sqrt(var2*ratio2^2)
        sd3<-sqrt(var3*ratio3^2)
        interval1<-c(ratio1-sd1*conv_tol, ratio1+sd1*conv_tol)
        interval2<-c(ratio2-sd2*conv_tol, ratio2+sd2*conv_tol)
        interval3<-c(ratio3-sd3*conv_tol, ratio3+sd3*conv_tol)
        overlap1<-overlap_check(interval1[1], interval1[2], interval2[1], interval2[2])
        overlap2<-overlap_check(interval1[1], interval1[2], interval3[1], interval3[2])
        overlap3<-overlap_check(interval2[1], interval2[2], interval3[1], interval3[2])
        cost_match<-c(overlap1, overlap2, overlap3)
        cost_conv<-all(cost_match)
        cost_archive[ql_iter-min_episode,]<-cost_match

        print(cost_match)

        if (ql_iter>max_episode) { break }
    }

    # record final search results
    solns<-list(soln1, soln2, soln3)
    ratios<-c(ratio1, ratio2, ratio3)
    vars<-c(var1, var2, var3)
    cost_archive<-cost_archive[complete.cases(cost_archive),]
    return(list(solns, ratios, vars, cost_archive))
}

## ql_episode subroutine represents a single episode of learning
# each wave explores a single path using one bundle of samples
# inputs: 
# r_map : list by indexer (one entry per node), within list, place list by adjacent nodes, place n x 3 matrix for BARratio, varBAR and varvarBAR for each bundle's evaluation
# q_map : list by indexer (one entry per node), each list entry  is vector of q values to neighbors

ql_episode<-function(q_map, r_map, epsilon) {

    # initialize this episode's search agent
    # agent holds info on current path and draws
    # at end, reprocess agent draws to determine var of each edge
    # report vars as rewards and update q map

    # create path and collect forward trajectory draws
    agent_path<-vector()
    agent_path[1]<-indexer_init
    iter_dummy<-0

    # create storage for forward draws and weights
    draws_F<-array(NA, dim=c(num_traj, 1))
    weights_F<-array(NA, dim=c(num_traj, 1))

    # initialize values 
    draws_F[,1]<-rnorm(num_traj, mean=mu0, sd=sig0)
    weights_F[,1]<-1/num_traj

    while(!identical(agent_path[length(agent_path)], indexer_target)) { #until the agent has reached the target, continue the search

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
        new_draws<-sapply(draws_F[,iter_dummy],smc_transition, indexer_next)
        draws_F<-cbind(draws_F, new_draws)

        # calculate particle weights
        p_weight<-get_p_weight(draws_F[,iter_dummy], weights_F[,iter_dummy], agent_path[iter_dummy], agent_path[iter_dummy+1])
        weights_F<-cbind(weights_F, p_weight)

        # continue particle propagation until path is complete
    }

    # propagate particles backwards
    draws_R<-array(NA, dim=dim(draws_F))
    weights_R<-array(NA, dim=dim(weights_F))
    draws_R[,ncol(draws_R)]<-rnorm(num_traj, mean=mu1, sd=sig1)
    weights_R[,ncol(weights_R)]<-1/num_traj

    for (i in (ncol(draws_R)-1):1) {
        draws_R[,i]<-sapply(draws_R[,i+1], smc_transition, agent_path[i])
        weights_R[,i]<-get_p_weight(draws_R[,i+1], weights_R[,i+1], agent_path[i+1], agent_path[i])
    }

    # process forward and reverse weighted trajectories to obtain ratio and var(ratio) estimates.
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

split_func<-function(x) { 
    if(any(is.na(x))) { out<-NA
    } else { out<-as.numeric(unlist(strsplit(x, split='_'))) }
    return(out)
}

overlap_check<-function(x1, x2, y1, y2) { return(x1<=y2 && y1<=x2) }

epsilon_update_cost<-function(epsilons, cost_match) {
    epsilon_delta<-0.005
    if (cost_match[1]) {
        epsilons[3]<-max(0, epsilons[3]-epsilon_delta)
        epsilons[1]<-min(1, epsilons[1]+epsilon_delta)
        epsilons[2]<-min(1, epsilons[2]+epsilon_delta)
    } else if (cost_match[2]) {
        epsilons[2]<-max(0, epsilons[2]-epsilon_delta)
        epsilons[1]<-min(1, epsilons[1]+epsilon_delta)
        epsilons[3]<-min(1, epsilons[3]+epsilon_delta)
    } else if (cost_match[3]) {
        epsilons[1]<-max(0, epsilons[1]-epsilon_delta)
        epsilons[3]<-min(1, epsilons[3]+epsilon_delta)
        epsilons[2]<-min(1, epsilons[2]+epsilon_delta)
    } else {
        epsilons[1]<-max(0, epsilons[1]-epsilon_delta)
        epsilons[2]<-max(0, epsilons[2]-epsilon_delta)
        epsilons[3]<-max(0, epsilons[3]-epsilon_delta)
        # epsilons[1]<-min(0, epsilons[1]+epsilon_delta)
        # epsilons[2]<-min(0, epsilons[2]+epsilon_delta)
        # epsilons[3]<-min(0, epsilons[3]+epsilon_delta)
    }
    return(epsilons)
}

# testing vars
# epsilon<-0
# alpha<-0.8
# gamma<-0.8