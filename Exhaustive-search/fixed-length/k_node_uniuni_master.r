# k_node_uniuni_master.r

# run path search on varBAR distance in a grid indexed by lambda and temperature
# initial and target distributions are normal distributions
# here we find paths of fixed length, k
# this allows us to adjust the path cost to account for the total number of edges in the path, which is a proxy of computational cost of traversing that path.

# solution based on http://stackoverflow.com/questions/1690347/shortest-path-with-a-fixed-number-of-edges

###############################
## PART 1: SET UP PARAMS ######
###############################

# set output directory
switch(Sys.info()[['sysname']],
    Windows = {setwd('D:/Dropbox/kshortest/')},
    Darwin = {setwd('/Users/rmuraglia/Dropbox/kshortest/')},
    Linux = {setwd('/dscrhome/rm204/kbest/')}
)

# load libraries silently
options(warn=-1)
suppressMessages(library(ggplot2))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
# suppressMessages(library(gplots))
options(warn=0)

# grab command line arguments, if there are any
args<-commandArgs(trailingOnly=TRUE)

if (length(args)==4) { # if proper number of params are provided, use them
    filID<-as.character(args[1])
    move.jump<-as.numeric(args[2])
    temp.numpoints<-as.numeric(args[3])
    lambda.numpoints<-as.numeric(args[4])
} else { # otherwise, use these defaults
    filID<-'knode03'
    move.jump<-4
    temp.numpoints<-10 
    lambda.numpoints<-6 
}

# set fixed params
reps<-5000
temp.min<-0.5
temp.max<-5
measure<-'asympvar'
nmax<-20000
dof<-c(TRUE, TRUE)

# end point distribution params
mu0<-0
mu1<-5
sig0<-1
sig1<-1

# set output path
filname<-paste('output/', filID, '-', move.jump, 'jump', sep='')

# set up grid, import methods
source('grid_import_uniuni.r')
source('dist_import_uniuni.r')
source('search_import_uniuni.r')

# collect draws
sample_list<-vector('list', length(indexer))
names(sample_list)<-indexer
for (i in 1:length(point.list)) {
    sample_list[[i]]<-collect_draws(point.list[[i]], reps)
}

# create adjacency lists
adj_list<-vector('list', length(indexer))
names(adj_list)<-indexer
for (i in 1:length(point.list)) {
    move.bool<-apply(move.universe, 1, is.valid.move, point.list[[i]])
    valid.moves<-move.universe[move.bool,]
    neighbors<-t(apply(valid.moves, 1, move.result, point.list[[i]]))
    adj_list[[i]]<-apply(neighbors, 1, paste, collapse='_')
}

##########
# PART 2: run the search for standard length 
##########

max_num_nodes<-lambda.numpoints + 2*(temp.numpoints-temp.ind)
dist_matrix<-array(NA, dim=c(length(indexer), length(indexer)))

search_out<-k_node_search(dist_matrix, max_num_nodes)

dist_matrix<-search_out[[1]]
unscaled_cost<-search_out[[2]]

########
# PART 3: visualize results
########

source('k_node_viz_uniuni.r')


