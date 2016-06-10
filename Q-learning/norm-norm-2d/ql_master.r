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
    num_episode<-200
    # min_episode<-50
}

# params for QL search
alpha<-0.8
gamma<-1
epsilon_init<-0
epsilon_tau<-170
# epsilon_delta<-0.0/min_episode # at the end of the minimum number of search episodes, set epsilon to 0.5
conv_tol<-0.0025

# set up grid
source('ql_grid_setup.r')

# load QL algorithms
source('ql_algorithm.r')

# num_traj<-200
# num_episode<-50
# epsilon_tau<-30

tim<-ql_full_history(q_map, r_map, alpha, gamma, epsilon_init, epsilon_tau)


library(ggplot2)
library(grid)
library(gridExtra)
library(gplots)

num_solns<-length(unique(tim[[1]]))

plot_df<-array(NA, dim=c(200, 3))
colnames(plot_df)<-c('iter', 'cost', 'soln')
plot_df[,1]<-c(1:200)
plot_df[,2]<-tim[[2]]
for (i in 1:num_solns) {
    inds<-which(sapply(tim[[1]], identical, unique(tim[[1]])[[i]]))
    plot_df[inds,3]<-i
}
plot_df<-as.data.frame(plot_df)
plot_df[,3]<-as.factor(plot_df[,3])

pcolors<-c('black', 'black', rich.colors(num_solns-2))
pcolors<-c('black',rich.colors(num_solns-1))

p0<-ggplot(data=plot_df, aes(x=iter, y=cost, colour=soln)) + geom_point() + theme(legend.position='top') + scale_colour_manual(values=pcolors)

p2<-ggplot(data=point_grid, aes(x=lambda, y=temperature)) + geom_path(data=as.data.frame(unique(tim[[1]])[[2]]), arrow=arrow(), colour=pcolors[2], size=3) + geom_point() + labs(title='soln 2')

p3<-ggplot(data=point_grid, aes(x=lambda, y=temperature)) + geom_path(data=as.data.frame(unique(tim[[1]])[[3]]), arrow=arrow(), colour=pcolors[3], size=3) + geom_point() + labs(title='soln 3')

p4<-ggplot(data=point_grid, aes(x=lambda, y=temperature)) + geom_path(data=as.data.frame(unique(tim[[1]])[[4]]), arrow=arrow(), colour=pcolors[4], size=3) + geom_point() + labs(title='soln 4')

p5<-ggplot(data=point_grid, aes(x=lambda, y=temperature)) + geom_path(data=as.data.frame(unique(tim[[1]])[[5]]), arrow=arrow(), colour=pcolors[5], size=3) + geom_point() + labs(title='soln 5')

p6<-ggplot(data=point_grid, aes(x=lambda, y=temperature)) + geom_path(data=as.data.frame(unique(tim[[1]])[[6]]), arrow=arrow(), colour=pcolors[6], size=3) + geom_point() + labs(title='soln 6')

p7<-ggplot(data=point_grid, aes(x=lambda, y=temperature)) + geom_path(data=as.data.frame(unique(tim[[1]])[[7]]), arrow=arrow(), colour=pcolors[7], size=3) + geom_point() + labs(title='soln 7')

p8<-ggplot(data=point_grid, aes(x=lambda, y=temperature)) + geom_path(data=as.data.frame(unique(tim[[1]])[[8]]), arrow=arrow(), colour=pcolors[8], size=3) + geom_point() + labs(title='soln 8')

p9<-ggplot(data=point_grid, aes(x=lambda, y=temperature)) + geom_path(data=as.data.frame(unique(tim[[1]])[[9]]), arrow=arrow(), colour=pcolors[9], size=3) + geom_point() + labs(title='soln 9')

p10<-ggplot(data=point_grid, aes(x=lambda, y=temperature)) + geom_path(data=as.data.frame(unique(tim[[1]])[[10]]), arrow=arrow(), colour=pcolors[10], size=3) + geom_point() + labs(title='soln 10')


outplot<-arrangeGrob(p0, p3, p4, p5, p6, p7, nrow=2, widths=c(1,1,1))

ggsave(file='maze-grid01-monitor-10.png', width=9, height=6, dpi=100, plot=outplot)
ggsave(file='maze-grid01-monitor-05.png', width=9, height=9, dpi=100, plot=outplot)
ggsave(file='maze-grid01-monitor-07.png', width=7, height=6, dpi=100, plot=outplot)

