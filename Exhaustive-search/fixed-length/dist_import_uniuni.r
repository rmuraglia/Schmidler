# dist_import_uniuni.r

#####################################
## PART 1: PARAMETER UPDATES ########
#####################################

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

#####################################
## PART 2: DISTANCE EVAL FUNC #######
#####################################

dist.move<-function(state1, state2) {
    if (measure=='asympvar') { dist<-asympvar_dist(state1, state2)
        } else { print('ERROR: invalid distance measure') }
    return(dist)
}

# params for asympvar_dist()

max_iter<-20 # max number of iterations for iterative ratio solver
tol<-1E-5 # tolerance for ratio convergence
init_r<-1 # some initial guess for the ratio

asympvar_dist<-function(state1, state2) {
    
    ## get a free energy estimate
    
    # grab samples from samplematrix
    state1key<-paste(state1,collapse='_')
    state1ind<-which(indexer==state1key)
    state2key<-paste(state2,collapse='_')
    state2ind<-which(indexer==state2key)
    draws1<-sample_list[[state1ind]]
    draws2<-sample_list[[state2ind]]

    # estimate ratio and free energy
    ratio_out<-ratio_estim(draws1, draws2, state1, state2)
    free_en<-log(ratio_out)

    ## calculate the ensemble average

    del_u1<-potential_func(draws1, state2[1], state2[2]) - potential_func(draws1, state1[1], state1[2])
    del_u2<-potential_func(draws2, state2[1], state2[2]) - potential_func(draws2, state1[1], state1[2])

    del_u<-c(del_u1, del_u2)

    # store other quantities
    M<-log(length(draws1)/length(draws2))
    N<-length(draws1)+length(draws2)

    ens_vals<-(2+2*cosh(free_en - del_u - M))^(-1)
    ens_mean<-mean(ens_vals)

    # get asymptotic variance
    asymp_out<-(1/N)*(ens_mean^(-1) - (N/length(draws2) + N/length(draws1)))

    # correct for number of samples to return intrinsic variance, unscaled by n
    dist_out<-asymp_out*N
    return(dist_out)
}

####################
# helper functions for asympvar_dist
####################

potential_func<-function(draw, lambda, temp) {
    mu_new<-mu_update(lambda, temp)
    sig_new<-sig_update(lambda, temp)
    potential<-(draw-mu_new)^2/(2*sig_new^2)
    return(potential)
}

unnorm_func<-function(draw, lambda, temp) {
    potential<-potential_func(draw, lambda, temp)
    unnorm<-exp(-potential)
    return(unnorm)
}

ratio_estim<-function(draws1, draws2, state1, state2) {
    
    ## precompute unnorm density ratios
    l1<-unnorm_func(draws1, state1[1], state1[2])/unnorm_func(draws1, state2[1], state2[2])
    l2<-unnorm_func(draws2, state1[1], state1[2])/unnorm_func(draws2, state2[1], state2[2])

    ## calculate sample balance variables
    s1<-length(draws1)/(length(draws1)+length(draws2))
    s2<-length(draws2)/(length(draws1)+length(draws2))

    ## estimate the ratio
    r_estim<-iterate_ratio(l1, l2, s1, s2)
    return(r_estim)
}

iterate_ratio<-function(l1, l2, s1, s2) {

    # initialize storage for ratio guess sequence
    r_vec<-rep(NA, length.out=max_iter)
    r_vec[1]<-init_r
    r_curr<-init_r
    r_prev<-r_curr+5
    iter<-1

    # keep iterating until we reach convergence or reach max iter
    while ((iter<max_iter) && abs(r_prev-r_curr)>tol) {
        r_prev<-r_vec[iter]
        r_curr<-update_ratio(r_prev, l1, l2, s1, s2)
        iter<-iter+1
        r_vec[iter]<-r_curr
    }

    # trim off unused storage and return final estimate
    r_vec<-r_vec[!is.na(r_vec)]
    return(r_vec[length(r_vec)])
}

update_ratio<-function(prev_r, l1, l2, s1, s2) {
    numerator<-(1/length(l2))*sum(l2/(s1*l2+s2*prev_r))
    denominator<-(1/length(l1))*sum(1/(s1*l1+s2*prev_r))
    new_r<-numerator/denominator
    return(new_r)
}
