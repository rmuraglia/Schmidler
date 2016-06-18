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
    Darwin = {setwd('/Users/rmuraglia/GitHub/Schmidler/Q-learning/norm-norm-2d/')},
    Linux = {setwd('/dscrhome/rm204/Q-learning/')}
)

print ('Reading parameters from command line if present...')
args<-commandArgs(trailingOnly=TRUE)
if (length(args)==1) { # if appropriate number, use
    repnum<-as.numeric(args[1])
    # move_jump<-as.numeric(args[1])
    # temp_numpoints<-as.numeric(args[2])
    # lambda_numpoints<-as.numeric(args[3])
    # temp_min<-as.numeric(args[4])
    # temp_max<-as.numeric(args[5])
    # temp_dof<-as.numeric(args[6])
    # num_traj<-as.numeric(args[7])
    # num_episode<-as.numeric(args[8])
} else { # else use following defaults
    move_jump<-1
    temp_numpoints<-5
    lambda_numpoints<-6
    temp_min<-0.5
    temp_max<-4
    temp_dof<-TRUE
    num_traj<-100
    num_episode<-600
}

move_jump<-1
temp_numpoints<-5
lambda_numpoints<-6
temp_min<-0.5
temp_max<-4
temp_dof<-TRUE
num_traj<-100
num_episode<-600

# end point distribution params
mu0<-0
mu1<-5
sig0<-1
sig1<-1

# params for sequential sampling 
nummetrostep<-25
metro_spread<-c(0.5)

# params for QL search
alpha<-0.8
gamma<-1
epsilon_init<-0
epsilon_tau<-floor(0.85*num_episode)
# epsilon_tau<-450 # determines episilon rate of change - end searching on epsilon_tau'th search episode

# set up grid
source('ql_grid_setup.r')

# load QL algorithms
source('ql_algorithm.r')

tim<-ql_full_history(q_map, r_map, alpha, gamma, epsilon_init, epsilon_tau)

source('ql_viz_fullhist.r')

# call ql_viz_fullshit.r


