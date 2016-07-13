# pcrooks_import.r

# import methods related to calculating the pcrooks distance. has depends on ising import for energy functions

get_p_weight<-function(x, prev_weight, indexer_prev, indexer_curr) {
    distr_prev<-as.numeric(unlist(strsplit(indexer_prev, split='_')))
    distr_curr<-as.numeric(unlist(strsplit(indexer_curr, split='_')))
    curr_log<-apply(x, 1, log_dens, distr_curr[1], distr_curr[2])
    prev_log<-apply(x, 1, log_dens, distr_prev[1], distr_prev[2])
    # curr_dens<-apply(x, 1, unnorm_dens, distr_curr[1], distr_curr[2])
    # prev_dens<-apply(x, 1, unnorm_dens, distr_prev[1], distr_prev[2])
    incr_weights<-exp(curr_log-prev_log)
    # incr_weights<-curr_dens/prev_dens
    unnorm_weights<-prev_weight*incr_weights
    p_weights<-unnorm_weights/sum(unnorm_weights)
    return(p_weights)
}

process_trajectories<-function(draws_F, draws_R, weights_F, weights_R, path) {

    # initialize storage for ratios and variances
    ratios<-rep(NA, length.out=length(path)-1)
    vars<-ratios

    for (i in 1:length(ratios)) {
        if (i==1) {
            curr_samp<-draws_F[ , ,i]
            next_samp<-draws_R[ , ,i+2]
        } else if (i==length(ratios)) {
            curr_samp<-draws_F[ , ,i-1]
            next_samp<-draws_R[ , ,i+1]
        } else {
            curr_samp<-draws_F[ , ,i-1]
            next_samp<-draws_R[ , ,i+2]
        }

        ratios[i]<-ratio_estim(curr_samp, next_samp, weights_F[,i], weights_R[,i+1], path[i], path[i+1])
        vars[i]<-var_estim(curr_samp, next_samp, weights_F[,i], weights_R[,i+1], path[i], path[i+1], ratios[i])
    }
    return(rbind(ratios, vars))
}

ratio_estim<-function(curr_samp, next_samp, curr_weight, next_weight, curr_indexer, next_indexer) {

    # NOTE: This calculates Z_next / Z_curr

    # get distribution info
    distr_curr<-as.numeric(unlist(strsplit(curr_indexer, split='_')))
    distr_next<-as.numeric(unlist(strsplit(next_indexer, split='_')))

    # precompute ratios of unnormalized densities
    # l_curr<-apply(curr_samp, 1, unnorm_dens, distr_next[1], distr_next[2])/apply(curr_samp, 1, unnorm_dens, distr_curr[1], distr_curr[2])
    # l_next<-apply(next_samp, 1, unnorm_dens, distr_next[1], distr_next[2])/apply(next_samp, 1, unnorm_dens, distr_curr[1], distr_curr[2])
    l_curr<-exp(apply(curr_samp, 1, log_dens, distr_next[1], distr_next[2]) - apply(curr_samp, 1, log_dens, distr_curr[1], distr_curr[2]))
    l_next<-exp(apply(next_samp, 1, log_dens, distr_next[1], distr_next[2]) - apply(next_samp, 1, log_dens, distr_curr[1], distr_curr[2]))

    # precompute sample balance
    s_curr<-length(l_curr)/(length(l_curr) + length(l_next))
    s_next<-length(l_next)/(length(l_curr) + length(l_next))

    # estimate ratio by iteration
    r_hat<-iterate_ratio(l_curr, l_next, s_curr, s_next, curr_weight, next_weight)
    return(r_hat)
}

iterate_ratio<-function(l_curr, l_next, s_curr, s_next, curr_weight, next_weight) {

    # iteration helpers
    max_iter<-20
    tol<-1E-5
    init_r<-1

    # initialize storage for ratio guess sequence
    r_vec<-rep(NA, length.out=max_iter)
    r_vec[1]<-init_r
    r_curr<-init_r
    r_prev<-r_curr+5
    iter<-1

    # keep iterating until we reach convergence or reach max iter
    while ((iter<max_iter) && abs(r_prev-r_curr)>tol) {
        
        # update best ratio guess
        r_prev<-r_vec[iter]

        # calculate new ratio
        numerator<-(1/length(l_curr)) * sum(curr_weight*l_curr/(s_next*l_curr+s_curr*r_prev))
        denominator<-(1/length(l_next)) * sum(next_weight/(s_next*l_next+s_curr*r_prev))
        r_curr<-numerator/denominator

        iter<-iter+1
        r_vec[iter]<-r_curr
    }

    # trim off unused storage and return final estimate
    r_vec<-r_vec[!is.na(r_vec)]
    return(r_vec[length(r_vec)])
}

var_estim<-function(curr_samp, next_samp, curr_weight, next_weight, curr_indexer, next_indexer, r_hat) {

    # get distn info
    distr_curr<-as.numeric(unlist(strsplit(curr_indexer, split='_')))
    distr_next<-as.numeric(unlist(strsplit(next_indexer, split='_')))

    # get delta us
    # del_u_curr<--log(apply(curr_samp, 1, unnorm_dens, distr_curr[1], distr_curr[2])) + log(apply(curr_samp,1 , unnorm_dens, distr_next[1], distr_next[2]))
    # del_u_next<--log(apply(next_samp, 1, unnorm_dens, distr_curr[1], distr_curr[2])) + log(apply(next_samp,1 , unnorm_dens, distr_next[1], distr_next[2]))
    del_u_curr<--apply(curr_samp, 1, log_dens, distr_curr[1], distr_curr[2]) + apply(curr_samp,1 , log_dens, distr_next[1], distr_next[2])
    del_u_next<--apply(next_samp, 1, log_dens, distr_curr[1], distr_curr[2]) + apply(next_samp,1 , log_dens, distr_next[1], distr_next[2])

    # get free energy and other quantities
    delta_f<--log(r_hat)
    M<-log(length(next_weight)/length(curr_weight))
    N<-length(curr_weight) + length(next_weight)

    # because weights collectively sum to 2, replace mean with sum/2
    ens_vals_curr<-curr_weight/(2+2*cosh(delta_f - del_u_curr - M))
    ens_vals_next<-next_weight/(2+2*cosh(delta_f - del_u_next - M))
    ens_mean<-sum(c(ens_vals_curr, ens_vals_next))/2

    # get variance
    var_out<-(1/N)*(ens_mean^(-1) - (N/length(curr_samp) + N/length(next_samp)))
    return(var_out)
}

