# yacas_precompute.r
# potentially not even a yacas style thing anymore? working with energies, can I just cirucmvent the need for symbolic math in R?

# args<-commandArgs(trailingOnly=TRUE)

# if (length(args)==3) {
#     num_distr<-as.numeric(args[1])
#     path_height<-as.numeric(args[2])
#     targ_dist<-as.character(args[3])
# }

# yacasID<-paste('yacas-import-', targ_dist, '-', path_height, 'height-', num_distr, 'distr', '.RData', sep='')

# set lambda and temperature values
lambda_min<-0
lambda_max<-1
lambda_points<-seq(from=lambda_min, to=lambda_max, length.out=num_distr)

init_temp<-1
T_get<-function(lambda, max_height) { round(sqrt(max_height^2 - max_height^2/0.5^2 * (lambda-0.5)^2) + init_temp, 3) }
temp_points<-T_get(lambda_points, path_height)

distr_mat<-rbind(lambda_points, temp_points)
distr_path<-as.list(as.data.frame(distr_mat))
distr_indexer<-sapply(distr_path, paste, collapse='_')

# yacas preliminaries
# library(Ryacas)
# Exp<-function(x) {exp(x)}
# x<-Sym('x')

# set initial distribution params
mu0<-0
sig0<-1
Z0<-sig0*sqrt(2*pi)
U0<-function(x) { 
    potential<-(x-mu0)^2/(2*sig0^2) 
    return(potential)
}

# set target distribution params
if (targ_dist=='norm') {
    mu1<-5
    sig1<-1
    U1<-function(x) {
        potential<-(x-mu1)^2/(2*sig1^2)
    }
} else if (targ_dist=='tdistn') {
    nu1<-1
    Z1<-(sqrt(nu1*pi)*gamma(nu1/2))/(gamma((nu1+1)/2))
    U1<-function(x) {
        potential<--log(dt(x, nu1)*Z1)
        return(potential)
    }
} else if (targ_dist=='bimod') {
    a<-0.5
    b<-14
    c<-64
    U1<-function(x) {
        potential<-(a*x^4 - b*x^2)/c
        return(potential)
    }
} else { 
    print('targ_dist flag should be one of "tdistn", "norm" or "bimod."')
    print('exiting R...')
    quit('no')
}

# set energy function
lambda_energy<-function(x, lambda, temp) {
    U_lam<-(1-lambda)*U0(x) + lambda*U1(x)
    return((1/temp)*U_lam)
}

# if bimodal target distribution, require samplepdf method for lambda=1 sampling
if (targ_dist=='bimod') {
    # library(Ryacas)
    # helper function
    endsign <- function(f, sign = 1) { #helper function to focus the search for uniroot
        b <- sign
        while (sign * f(b) < 0) { b <- b + sign*5 } #can tune how big to make search, depending on how spread the distribution is
        #while (sign * f(b) < 0) b <- 10 * b #original search expansion
        return(b)
    }

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

    bimod_unnorm<-function(x) { exp(-lambda_energy(x, lambda_max, init_temp)) }

    bimod_Z<-integrate(bimod_unnorm, -Inf, Inf)$value

    bimodpdf<-function(x) { bimod_unnorm(x)/bimod_Z }
}

# save.image(file=yacasID)