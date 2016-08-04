# precompute_draws_unibi.r

##################
# independent sampler engine
##################

# helper function
endsign <- function(f, sign = 1) { #helper function to focus the search for uniroot
    b <- sign
    while (sign * f(b) < 0) { b <- b + sign*5 } #can tune how big to make search, depending on how spread the distribution is
    #while (sign * f(b) < 0) b <- 10 * b #original search expansion
    return(b)
}

# get normalization constant
true_Z<-function(lambda, temperature) {

    # get unnorm expression based on indexer lookup
    key<-paste(lambda, temperature, sep='_')
    ind<-which(indexer==key)
    expr<-unnorm_expr_vec[ind]

    # parse function for calculation
    qfunc<-function(x){eval(parse(text=expr))}
    Z<-integrate(qfunc, -Inf, Inf)
    return(Z$value)
}

# sample from any pdf
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

################
# collect the samples
################

sample_matrix<-array(NA, dim=c(num_draws, length(indexer)))
colnames(sample_matrix)<-indexer

print(paste('there are a total of ', length(indexer), ' distributions to be sampled', sep=''))
print(paste('each distribution will be sampled ', num_draws, ' times', sep=''))

# set up checkpoints for run-time diagonistics
check01<-ceiling(length(indexer)*0.01)
check05<-ceiling(length(indexer)*0.05)
check10<-ceiling(length(indexer)*0.1)
checkcount<-0
checkprint<-c(1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
checkind<-1

for (i in 1:length(indexer)) {

    # get draws
    tmp_draws<-rep(NA, length.out=num_draws)
    Z<-true_Z(point_list[[i]][1], point_list[[i]][2])
    drawfunc<-function(x){eval(parse(text=unnorm_expr_vec[i]))/Z}
    tmp_draws<-samplepdf(num_draws, drawfunc)
    sample_matrix[,i]<-tmp_draws

    # print progress
    checkcount<-checkcount+1
    if (checkcount<=check05 && checkcount%%check01==0) {
        print(date())
        print(paste(checkprint[checkind], '% done sampling', sep=''))
        checkind<-checkind+1
    }
    if (checkcount%%check10==0) {
        print(date())
        print(paste(checkprint[checkind], '% done sampling', sep=''))
        checkind<-checkind+1
    }
}

write.table(sample_matrix, file=paste('samplematrix-', gridID, '-', num_draws, 'numdraws', '.txt', sep=''), append=F, quote=F, sep='\t', col.names=T, row.names=F)