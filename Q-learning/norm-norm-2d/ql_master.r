# ql_master.r

# master file controlling Q-learning search
# initial distribution: N(0,1)
# target distribution: N(5,1)
# add additional info here: distance metric, sampling method

#############
# PART 1: set up and initilization
#############

print('Begin R script')
print(date())

# set working directory
switch(Sys.info()[['sysname']],
    Windows = {setwd('D:/work-duke/Schmidler/Q-learning/norm-norm-2d')},
    Darwin = {setwd('/Users/rmuraglia/Dropbox/pathsearch-exhaustive/uni-bi-symmetric/')},
    Linux = {setwd('/dscrhome/rm204/Q-learning/')}
)

print ('Reading parameters from command line if present...')
args<-commandArgs(trailingOnly=TRUE)
if (length(args)==7) { # if appropriate number, use
    move_jump<-as.numeric(args[1])
    temp_numpoints<-as.numeric(args[2])
    lambda_numpoints<-as.numeric(args[3])
    temp_min<-as.numeric(args[4])
    temp_max<-as.numeric(args[5])
    temp_dof<-as.numeric(args[6])
    num_traj<-as.numeric(args[7])
} else { # else use following defaults
    move_jump<-1
    temp_numpoints<-5
    lambda_numpoints<-6
    temp_min<-0.5
    temp_max<-4
    temp_dof<-TRUE
    num_traj<-500
}

# end point distribution params
mu0<-0
mu1<-5
sig0<-1
sig1<-1

# set file ID


# set up grid