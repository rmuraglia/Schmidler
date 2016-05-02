# ql_grid_setup.r

# set up grid for Q-learning search in lambda-T parameter space
# to be called by ql_master.r - has depends on params defined in that file

##################
# set up grid
##################

lambda_min<-0
lambda_max<-1

lambda_points<-seq(from=lambda_min, to=lambda_max, length.out=lambda_numpoints)
temp_points<-seq(from=temp_min, to=temp_max, length.out=temp_numpoints) # note: temp means temperature, not temporary

# list of length 2 containing param values
point_lists<-list(lambda_points, temp_points)
names(point_lists)<-c('lambda', 'temperature')

# matrix containing combinations of param values
point_grid<-expand.grid(point_lists)

# list of combinations of param values
point_list<-as.list(as.data.frame(t(point_grid)))

# unique string identifiers for each param combination
indexer<-apply(point_grid, 1, paste, collapse='_')

################
# define intermediate distributions
################

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

# make a function that will take lambda, T and a draw to generate a potential energy, or an unnormalized density