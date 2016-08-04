# yacas_precompute.r

# get grid params from cmd line
args<-commandArgs(trailingOnly=TRUE)

if (length(args)==4) {
    temp_numpoints<-as.numeric(args[1])
    lambda_numpoints<-as.numeric(args[2])
    temp_min<-as.numeric(args[3])
    temp_max<-as.numeric(args[4])
} else {
    print('ERROR: must specify grid parameters')
    quit('no')
}

yacasID<-paste('yacas-import-', temp_numpoints,'tpts-', lambda_numpoints, 'lpts-', temp_min, 'tmin-', temp_max, 'tmax', '.RData', sep='')

## set fixed distn params
# unimodal distn params
mu0<-0
sig0<-1

#bimodal distn params
a<-0.5
b<-14
c<-64

## set up point grid
lambda_min<-0
lambda_max<-1
lambda_points<-seq(from=lambda_min, to=lambda_max, length.out=lambda_numpoints)
temp_points<-seq(from=temp_min, to=temp_max, length.out=temp_numpoints)

point_lists<-list(lambda_points, temp_points)
names(point_lists)<-c('lambda', 'temperature')
point_grid<-expand.grid(point_lists)
point_list<-as.list(as.data.frame(t(point_grid)))
indexer<-apply(point_grid, 1, paste, collapse='_')

## get unnormalized density expressions
library(Ryacas)
Exp<-function(x) {exp(x)}
x<-Sym('x')

U0<-function(mu, sigma) {
    potential<-(x-mu)^2/(2*sigma^2)
    return(potential)
}

U1<-function(A, B, C) {
    potential<-(A*x^4 - B*x^2)/C
    return(potential)
}

Umaker<-function(lambda) {
    newU<-(1-lambda)*U0(mu0, sig0) + lambda*U1(a,b,c)
    return(newU)
}

Qmaker<-function(lambda, temperature) {
    newQ<-exp(-(1/temperature)*Umaker(lambda))
    return(newQ)
}

unnorm_expr_vec<-rep(NA, length.out=length(point_list))
for (i in 1:length(point_list)) {
    unnorm_expr_vec[i]<-Qmaker(point_list[[i]][1], point_list[[i]][2])
}

# save to R object
save.image(file=yacasID)
