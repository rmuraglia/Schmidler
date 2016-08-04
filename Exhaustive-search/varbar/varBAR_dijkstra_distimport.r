# varBAR_dijkstra_distimport.r

# Define methods for distance evaluation for the varBAR tests

# library(Ryacas)

Qget<-function(potential,temperature) {
    integrand<-exp(-(1/temperature)*potential)
    Q<-integrate(function(x){eval(parse(text=integrand[1]))},-Inf,Inf)
    return(Q)
}

pdfmaker<-function(lambda, temperature) {
    U<-Umaker(lambda)
    Q<-Qget(U, temperature)
    newpdf<-exp(-(1/temperature)*U)/Q$value
    return(newpdf)
}

dist.move<-function(state1, state2) {
    if (measure=='empvar') { dist<-empvar_dist(state1, state2)
    } else if (measure=='asympvar') { dist<-asympvar_dist(state1, state2)
    } else { print('ERROR: invalide distance meaure') }
    return(dist)
}

empvar_dist<-function(state1, state2) {
    estim_vec<-rep(NA, length.out=numbootstrap)

    # set state info
    state1key<-paste(state1,collapse='_')
    state1ind<-which(indexer==state1key)
    state2key<-paste(state2,collapse='_')
    state2ind<-which(indexer==state2key)
    
    #define unnormalized distributions for states
    lam0<-state1[1]; beta0<-1/state1[2]
    lam1<-state2[1]; beta1<-1/state2[2]
    q0<-exp(-beta0*Umaker(lam0))
    q1<-exp(-beta1*Umaker(lam1))
    q0func<-function(x){eval(parse(text=q0[1]))}
    q1func<-function(x){eval(parse(text=q1[1]))}

    estim_vec<-replicate(numbootstrap, bootfunc(samplematrix, state1ind, state2ind, q0func, q1func))
    empvar_out<-var(estim_vec)/mean(estim_vec)^2
    return(empvar_out)
}

bootfunc<-function(samplematrix, state1ind, state2ind, q0func, q1func) {
    bootrand<-sample(1:reps, size=reps, replace=T)
    bootdraws1<-samplematrix[bootrand, state1ind]
    bootdraws2<-samplematrix[bootrand, state2ind]
    bootestim<-ratio_estim(q0func, q1func, bootdraws1, bootdraws2)
    return(bootestim)
}

max_iter<-20
tol<-1E-5
init_r<-1

ratio_estim<-function(func0, func1, draws1, draws2) {
    # precompute ratios
    l1<-func0(draws1)/func1(draws1)
    l2<-func0(draws2)/func1(draws2)
    s1<-length(draws1)/(length(draws1)+length(draws2))
    s2<-length(draws2)/(length(draws1)+length(draws2))

    r_estim<-iterate_ratio(max_iter, tol, init_r, l1, l2, s1, s2)
    return(r_estim)
}

iterate_ratio<-function(max_iter, tol, init_r, l1, l2, s1, s2) {
    r_vec<-rep(NA, length.out=max_iter)
    r_vec[1]<-init_r
    r_curr<-init_r
    r_prev<-r_curr+5
    iter<-1

    while ((iter<max_iter) && abs(r_prev-r_curr)>tol) {
        r_prev<-r_vec[iter]
        r_curr<-update_ratio(r_prev, l1, l2, s1, s2)
        iter<-iter+1
        r_vec[iter]<-r_curr
    }
    r_vec<-r_vec[!is.na(r_vec)]
    return(r_vec[length(r_vec)])
}

update_ratio<-function(prev_r, l1, l2, s1, s2) {
    numerator<-(1/length(l2))*sum(l2/(s1*l2+s2*prev_r))
    denominator<-(1/length(l1))*sum(1/(s1*l1+s2*prev_r))
    new_r<-numerator/denominator
    return(new_r)
}

asympvar_dist<-function(state1, state2) {
    ## get free en estim
    # set state info
    state1key<-paste(state1,collapse='_')
    state1ind<-which(indexer==state1key)
    state2key<-paste(state2,collapse='_')
    state2ind<-which(indexer==state2key)
    draws1<-samplematrix[,state1ind]
    draws2<-samplematrix[,state2ind]
    
    #define unnormalized distributions for states
    lam0<-state1[1]; beta0<-1/state1[2]
    lam1<-state2[1]; beta1<-1/state2[2]
    q0<-exp(-beta0*Umaker(lam0))
    q1<-exp(-beta1*Umaker(lam1))
    q0func<-function(x){eval(parse(text=q0[1]))}
    q1func<-function(x){eval(parse(text=q1[1]))}

    # ratio estimation
    ratio_out<-ratio_estim(q0func, q1func, draws1, draws2)
    free_en<-log(ratio_out)

    # calculate that ensemble average
    u0<-beta0*Umaker(lam0)
    u1<-beta1*Umaker(lam1)
    u0func<-function(x){eval(parse(text=u0[1]))}
    u1func<-function(x){eval(parse(text=u1[1]))}

    del_u<-c(u1func(draws1)-u0func(draws1), u1func(draws2)-u0func(draws2))
    M<-log(length(draws1)/length(draws2))
    N<-length(draws1)+length(draws2)

    ens_vals<-(2+2*cosh(free_en-del_u-M))^-1
    ens_mean<-mean(ens_vals)

    asymp_out<-(1/N)*(ens_mean^-1 - (N/length(draws2)+N/length(draws1)))
    return(asymp_out)
}

##
# legacy stuff for viz
##

endsign <- function(f, sign = 1) { #helper function to focus the search for uniroot
    b <- sign
    while (sign * f(b) < 0) { b <- b + sign*5 } #can tune how big to make search, depending on how spread the distribution is
    #while (sign * f(b) < 0) b <- 10 * b #original search expansion
    return(b)
}
#initially focus search in [-1,1], but if value isn't contained within, expand by adding +/- 5

samplepdf <- function(n, pdf, ..., spdf.lower = -Inf, spdf.upper = Inf) {
    vpdf <- function(v) sapply(v, pdf, ...)  # vectorize
    cdf <- function(z) integrate(vpdf, spdf.lower, z)$value
    invcdf <- function(u) {
        subcdf <- function(t) cdf(t) - u
        if (spdf.lower == -Inf)
        spdf.lower <- endsign(subcdf, -1)
        if (spdf.upper == Inf)
        spdf.upper <- endsign(subcdf)
        return(uniroot(subcdf, c(spdf.lower, spdf.upper))$root)
    }
    sapply(runif(n,min=0.01,max=0.99), invcdf)
}

ratio.estim<-function(state1, state2, draws) {  
    #get indexer values for states
    state1key<-paste(state1,collapse='_')
    state1ind<-which(indexer==state1key)
    state2key<-paste(state2,collapse='_')
    state2ind<-which(indexer==state2key)
    
    #define unnormalized distributions for states
    lam0<-state1[1]; beta0<-1/state1[2]
    lam1<-state2[1]; beta1<-1/state2[2]
    q0<-exp(-beta0*Umaker(lam0))
    q1<-exp(-beta1*Umaker(lam1))
    q0func<-function(x){eval(parse(text=q0[1]))}
    q1func<-function(x){eval(parse(text=q1[1]))}
    
    #get BAR estimate

    # get draws from states
    draws1<-draws[,state1ind]
    draws2<-draws[,state2ind]

    #precompute ratios of densities and size params
    l1<-q0func(draws1)/q1func(draws1)
    l2<-q0func(draws2)/q1func(draws2)
    s1<-length(draws1)/(length(draws1)+length(draws2))
    s2<-length(draws2)/(length(draws1)+length(draws2))
    
    #estimate ratio iteratively
    r.estim<-iterate_ratio(max_iter, tol, init_r, l1, l2, s1, s2)
    return(r.estim)
}
