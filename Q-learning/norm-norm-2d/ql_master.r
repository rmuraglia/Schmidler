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
if (length(args)==7) { # if appropriate number, use
    move_jump<-as.numeric(args[1])
    temp_numpoints<-as.numeric(args[2])
    lambda_numpoints<-as.numeric(args[3])
    temp_min<-as.numeric(args[4])
    temp_max<-as.numeric(args[5])
    temp_dof<-as.numeric(args[6])
    min_episode<-as.numeric(args[7])
} else { # else use following defaults
    move_jump<-1
    temp_numpoints<-5
    lambda_numpoints<-6
    temp_min<-0.5
    temp_max<-4
    temp_dof<-TRUE
    min_episode<-30
}

# params for QL search
alpha<-0.8
gamma<-1
epsilon_init<-0
epsilon_delta<-0.0/min_episode # at the end of the minimum number of search episodes, set epsilon to 0.5
conv_tol<-0.0025

# set up grid
source('ql_grid02_setup.r')

# load QL algorithms
source('ql_algorithm.r')

# num_traj<-200
# num_episode<-50
# epsilon_tau<-30

tim<-ql_three_chain(q_map, r_map, alpha, gamma, min_episode, conv_tol)

library(ggplot2)
library(grid)

split_func<-function(x) { 
    if(any(is.na(x))) { out<-NA
    } else { out<-as.numeric(unlist(strsplit(x, split='_'))) }
    return(out)
}

path_df<-as.data.frame(t(sapply(tim[[1]], split_func)))
colnames(path_df)<-c('lambda', 'temperature')

p0<-ggplot(data=point_grid, aes(x=lambda, y=temperature)) + geom_point() + geom_path(data=path_df, arrow=arrow()) + labs(title='converged in 30+55 steps')

# p0<-ggplot(data=point_grid, aes(x=lambda, y=temperature)) + geom_point() + geom_path(data=path_df, arrow=arrow()) + labs(title='failed to converge after 400 steps')

ggsave(file='maze-grid02-3chains-05.png', width=6, height=6, dpi=100, plot=p0)

# tim<-ql_search(q_map, r_map, alpha, gamma, epsilon_init, min_episode, conv_tol)


# soln<-ql_path_soln(tim[[1]])
# soln
# ql_path_cost(tim[[2]], soln)


# ql_path_cost2(tim[[2]], soln)