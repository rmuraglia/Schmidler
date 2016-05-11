# ql_algorithm.r

# to be called by ql_master.r
# define Q-learning algorithm. adapted to minimize cost instead of maximize reward

ql_search<-function(q_map, r_map, alpha, gamma, epsilon_init, epsilon_tau) {

    # carry out a preset number of search episodes
    delta_epsilon<-1/epsilon_tau
    maps<-list(q_map, r_map)

    for (i in 1:epsilon_tau) {
        print(i)
        maps<-ql_episode(maps[[1]], maps[[2]], epsilon_init-(i-1)*delta_epsilon)
    }
    for (i in (epsilon_tau+1):num_episode) {
        print(i)
        maps<-ql_episode(maps[[1]], maps[[2]], 1)
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
    agent_draws<-array(NA, dim=c(num_traj, 1))
    agent_draws[,1]<-rnorm(num_traj, mean=mu0, sd=sig0)
    agent_incr_weights<-array(NA, dim=c(num_traj, 1))
    agent_norm_weights<-array(1/num_traj, dim=c(num_traj, 1))

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

        # propagate particles
        new_draws<-sapply(agent_draws[,iter_dummy],sis_transition, indexer_next)
        agent_draws<-cbind(agent_draws, new_draws)

        # calculate weights for ratio estimation
        agent_incr_weights<-cbind(agent_incr_weights, get_incr_weight(agent_draws[,iter_dummy], agent_path[iter_dummy], agent_path[iter_dummy+1]))
        agent_norm_weights<-cbind(agent_norm_weights, (agent_norm_weights[,iter_dummy]*agent_incr_weights[,iter_dummy+1])/sum(agent_norm_weights[,iter_dummy]*agent_incr_weights[,iter_dummy+1]))

        # bootstrap to get variance of ratio estimation
        estim_vec<-replicate(numbootstrap, boot_var_ratio(agent_norm_weights[,iter_dummy], agent_incr_weights[,iter_dummy+1]))
        action_reward<-var(log(estim_vec))
        reward_entry<-c(action_reward, 1) # for SIS case, start with uninformative 1 variance of variance - all estimates of reward are considered equally important

        # update r_map
        r_map[[curr_index]][[next_index]]<-rbind(r_map[[curr_index]][[next_index]], reward_entry)

        # update q_map
        q_entry<-as.numeric(q_map[[curr_index]][next_index, 2])
        r_entry<-mean(r_map[[curr_index]][[next_index]][,1] , na.rm=TRUE)
        q_prime_ind<-which(indexer==indexer_next)
        # q_prime_val<-min(as.numeric(q_map[[q_prime_ind]][,2]))
        # q_map[[curr_index]][next_index, 2]<-q_score(q_entry, r_entry, q_prime_val)
        q_prime_val<-max(as.numeric(q_map[[q_prime_ind]][,2]))
        q_map[[curr_index]][next_index, 2]<-q_score(q_entry, -r_entry, q_prime_val)

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

# testing vars
# epsilon<-0
# alpha<-0.8
# gamma<-0.8