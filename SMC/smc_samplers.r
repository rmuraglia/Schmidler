# smc_samplers.r

# contains methods for nonequilibrium sampling (with and without resampling)
# along with methods to convert (weighted) noneq draws to ratio estimates by AIS, sBAR, Crooks and pCrooks

#############
# PART 1: trajectory generation
#############

collect_noneq_traj<-function(prc_F) {
    if (prc_F==0) { # limiting case 1: no forward trajectories
        traj_F<-list(NA, NA)
        traj_R<-traj_gen_R(numtraj, distr_path)
    } else if (prc_F==1) { # other limiting case: no reverse trajectories
        traj_F<-traj_gen_F(numtraj, distr_path)
        traj_R<-list(NA, NA)
    } else {
        traj_F<-traj_gen_F(ceiling(numtraj*prc_F), distr_path)
        traj_R<-traj_gen_R(ceiling(numtraj*(1-prc_F)), distr_path)
    }
    return(list(traj_F, traj_R))
}

traj_gen_F<-function(bundle_size, distr_path) {

    # initialize storage
    samples_out<-array(NA, dim=c(bundle_size, length(distr_path)))
    incr_weights<-array(NA, dim=c(bundle_size, length(distr_path)-1))

    # get initial samples - in all cases here, forward draws originate from standard normal
    samples_out[,1]<-rnorm(bundle_size, mean=mu0, sd=sig0)

    # for each transition, place particle in new distn and mix
    for (i in 2:(length(distr_path))) {

        # get the density object
        curr_energy<-function(x) { lambda_energy(x, distr_path[[i-1]][1], distr_path[[i-1]][2]) }
        target_energy<-function(x) { lambda_energy(x, distr_path[[i]][1], distr_path[[i]][2]) }

        # get incremental weight of that transition
        # note: depends only on pre-mixed sample
        incr_weights[,i-1]<-get_incr_weight(samples_out[,i-1], target_energy, curr_energy)

        if (smc_meth=='seqbar') { # if it seqBAR, then resample

            # multinomial resample based on weights
            resample<-sample(samples_out[,i-1], size=bundle_size, replace=T, prob=incr_weights[,i-1])

            # propagate (mix to 'next' distn)
            samples_out[,i]<-sapply(resample, transition_func, target_energy)
        } else { # other if it is AIS, Crooks or pCrooks, skip resampling
            samples_out[,i]<-sapply(samples_out[,i-1], transition_func, target_energy)
        }
    } # close FOR

    # return draws and weights
    return(list(samples_out, incr_weights))
}

traj_gen_R<-function(bundle_size, distr_path) {

    # initalize storage
    samples_out<-array(NA, dim=c(bundle_size, length(distr_path)))
    incr_weights<-array(NA, dim=c(bundle_size, length(distr_path)-1))

    # get initial samples - depends on test case
    if (targ_dist=='tdistn') {
        samples_out[,length(distr_path)]<-rt(bundle_size, df=nu1)
    } else if (targ_dist=='norm') {
        samples_out[,length(distr_path)]<-rnorm(bundle_size, mean=mu1, sd=sig1)
    } else if (targ_dist=='bimod') {
        samples_out[,length(distr_path)]<-samplepdf(bundle_size, bimodpdf)
    } else { 
        print('targ_dist flag should be one of "tdistn", "norm" or "bimod."')
        print('exiting R...')
        quit('no')
    }

    # for each transition, place particle in new distn and mix
    for (i in (length(distr_path)-1):1) {

        curr_energy<-function(x) { lambda_energy(x, distr_path[[i+1]][1], distr_path[[i+1]][2]) }
        target_energy<-function(x) { lambda_energy(x, distr_path[[i]][1], distr_path[[i]][2]) }

        incr_weights[,i]<-get_incr_weight(samples_out[,i+1], target_energy, curr_energy)

        if (smc_meth=='seqbar') {
            resample<-sample(samples_out[,i+1], size=bundle_size, replace=T, prob=incr_weights[,i])
            samples_out[,i]<-sapply(resample, transition_func, target_energy)
        } else {
            samples_out[,i]<-sapply(samples_out[,i+1], transition_func, target_energy)
        }
    }
    return(list(samples_out, incr_weights))
}

## helper functions

get_incr_weight<-function(draws, targ_en, sampl_en) {
    incr_weights<-exp(-targ_en(draws))/exp(-sampl_en(draws))
}

# metropolis transition
metrospread<-0.5

transition_func<-function(curr_draw, energy_func) {
    A<-curr_draw
    for (i in 1:nummetrostep) {
        for (j in 1:length(metrospread)){ # allow for possibility of sequence of metropolis steps
            B<-metro_func(A, energy_func, metrospread[j])
            A<-B
        }
    }
    return(B)
}

metro_func<-function(curr_draw, energy_func, spread) {
    trial_draw<-rnorm(1, mean=curr_draw, sd=spread)
    
    curr_energy<-energy_func(curr_draw)
    trial_energy<-energy_func(trial_draw)
    energy_diff<-curr_energy - trial_energy
    rand_uni<-log(runif(n=1, min=0, max=1))
    if (rand_uni <= energy_diff) { return(trial_draw) 
    } else { return(curr_draw) }
}

###################
# PART 2: ratio estimation
###################

smc_estim<-function(draws_F, draws_R, weights_F, weights_R) {
    if (smc_meth=='ais') {
        if (prc_F==1) {
            estim_out<-ais_estim(weights_F)
        } else if (prc_F==0) {
            estim_out<-ais_estim(weights_R)
        } else {
            estim_F<-ais_estim(weights_F)
            estim_R<-ais_estim(weights_R)
            estim_C<-(estim_F[1]*(1/estim_F[2]) + estim_R[1]*(1/estim_R[2]))/(1/estim_F[2] + 1/estim_R[2])
            estim_out<-c(estim_C, (estim_F[2]+estim_R[2])/2)
        }
    } else if (smc_meth=='seqbar') {
        # calculate one ratio for each transition - per-step work
        ratios_out<-rep(NA, length.out=length(distr_path)-1)
        for (i in 1:length(ratios_out)) {
            if (prc_F==1) {
                work_F<--log(weights_F[,i])
                work_R<-get_sbar_alt_work(draws_F, i, 'F')
            } else if (prc_F==0) {
                work_F<-get_sbar_alt_work(draws_R, i, 'R') 
                work_R<--log(weights_R[,i])
            } else {
                work_F<-c(-log(weights_F[,i]), get_sbar_alt_work(draws_R, i, 'R')) 
                work_R<-c(-log(weights_R[,i]), get_sbar_alt_work(draws_F, i, 'F'))
            }
            ratios_out[i]<-bridge_estim(work_F, work_R, 1, 1)
        }
        estim_out<-prod(ratios_out)
    } else if (smc_meth=='crooks') {
        work_F<--log(weights_F)
        work_R<--log(weights_R)
        tot_work_f<-apply(work_F, 1, sum)
        tot_work_r<-apply(work_R, 1, sum)
        estim_out<-bridge_estim(tot_work_f, tot_work_r, 1, 1)
    } else if (smc_meth=='pcrooks') {
        # get normalized weights
        norm_F<-get_norm_weight(weights_F, 'F')
        norm_R<-get_norm_weight(weights_R, 'R')
        # get one ratio for each transition
        ratios_out<-rep(NA, length.out=length(distr_path)-1)
        for (i in 1:length(ratios_out)) {
            if (i==1) {
                work_F<--log(weights_F[,i])
                work_R<--log(weights_R[,i+1])
            } else if (i==length(ratios_out)) {
                work_F<--log(weights_F[,i-1])
                work_R<--log(weights_R[,i])
            } else { 
                work_F<--log(weights_F[,i-1])
                work_R<--log(weights_R[,i+1])
            }
            ratios_out[i]<-bridge_estim(work_F, work_R, norm_F[,i], norm_R[,i])
        }
        estim_out<-prod(ratios_out)
    } else {
        print('smc_meth should be one of "ais", "seqbar", "crooks" or "pcrooks."')
        print('exiting R...')
        quit('no')
    }
    return(estim_out)
}

ais_estim<-function(weights) {
    traj_weight<-apply(weights, 1, prod)
    r<-mean(traj_weight)
    v<-var(traj_weight)
    return(c(r, v))
}

get_sbar_alt_work<-function(draws, ind, FR) {

    en_0<-function(x) { lambda_energy(x, distr_path[[ind]][1], distr_path[[ind]][2]) }
    en_1<-function(x) { lambda_energy(x, distr_path[[ind+1]][1], distr_path[[ind+1]][2]) }

    if (FR=='F') {
        weights<-get_incr_weight(draws[,ind+1], en_0, en_1)
    } else if (FR=='R') {
        weights<-get_incr_weight(draws[,ind], en_1, en_0)
    } else{
        print('Something went wrong. This function should not have been called')
        print('FR flag should be "F" or "R"')
    }
    return(-log(weights))
}

get_norm_weight<-function(weights, FR) {
    norm_weights<-array(NA, dim=dim(weights))
    if (FR=='F') {
        norm_weights[,1]<-1/nrow(weights)
        for (i in 1:(ncol(weights)-1)) {
            norm_weights[,i+1]<-(norm_weights[,i]*weights[,i])/sum(norm_weights[,i]*weights[,i])
        }
    } else if (FR=='R') {
        norm_weights[,ncol(weights)]<-1/nrow(weights)
        for (i in (ncol(weights)-1):1) {
            norm_weights[,i]<-(norm_weights[,i+1]*weights[,i+1])/sum(norm_weights[,i+1]*weights[,i+1])
        }
    } else {
        print('Something went wrong. This function should not have been called')
        print('FR flag should be "F" or "R"')
    }
    return(norm_weights)
}

# ratio estimation params
max_iter<-20
tol<-1E-5
init_r<-1

bridge_estim<-function(work_f, work_r, norm_f, norm_r) {

    # initialize storage for ratio guess sequence
    r_vec<-rep(NA, length.out=max_iter)
    r_vec[1]<-init_r
    r_curr<-init_r
    r_prev<-r_curr+5 # this can be anything - just set r_curr != r_prev
    iter<-1

    # keep iterating until we reach convergence or reach max iter
    while ((iter<max_iter) && abs(r_prev-r_curr)>tol) {
        r_prev<-r_vec[iter]

        numerator<-(1/length(work_f))*sum(norm_f/(1+exp(work_f)*(length(work_f)/length(work_r))*r_prev))
        denominator<-(1/length(work_r))*sum(norm_r*exp(-work_r)/(1+exp(-work_r)*(length(work_f)/length(work_r))*r_prev))
        r_curr<-numerator/denominator
        iter<-iter+1
        r_vec[iter]<-r_curr
    }

    # trim off unused storage and return final estimate
    r_vec<-r_vec[!is.na(r_vec)]
    return(r_vec[length(r_vec)])
}











