# varBAR_dijkstra_uniuni_master.r

# run dijkstra search for varBAR distance in a grid indexed by lambda and temperature
# initial and target distributions are normal distributions

###############################
## PART 1: SET UP PARAMS ######
###############################

# set output directory
switch(Sys.info()[['sysname']],
    Windows = {setwd('D:/Dropbox/BAR estimator/')},
    Darwin = {setwd('/Users/ryan/Dropbox/BAR estimator')}
)

# load libraries silently
options(warn=-1)
suppressMessages(library(ggplot2))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(Ryacas))
suppressMessages(library(gplots))
options(warn=0)

# grab command line arguments, if there are any
args<-commandArgs(trailingOnly=TRUE)

if (length(args)==2) { #if the proper number of params are provided, use them
    filID<-as.character(args[1])
    move.jump<-as.numeric(args[2])
} else { #otherwise, use the ones below
    filID<-'varBARasymp-uni-uni-01'
    move.jump<-5
}

# set fixed params
reps<-1000
temp.min<-0.5
temp.max<-5
temp.numpoints<-10
lambda.numpoints<-6
# numbootstrap<-200 # only used for empvar by bootstrap
measure<-'asympvar'

# set filename prefix
filname<-paste('asympvar/', filID, 'b-', move.jump, 'jump', sep='')

################
## PART 2: Optimal path in full grid
################

source('varBAR_dijkstra_distimport.r')

# uni-uni target dist defs

Exp<-function(x) {exp(x)}
x<-Sym('x')

mu0<-0
mu1<-5

U<-function(mu) {
    potential<-(x-mu)^2/2
    return(potential)
}

Umaker<-function(lambda) {
    newU<-(1-lambda)*U(mu0) + lambda*U(mu1)
    return(newU)
}

# end target dist defs

dof<-c(TRUE, TRUE)
source('varBAR_dijkstra_import.r')
source('varBAR_dijkstra_precomputedraws.r')

samplematrix<-collect.draws(indexer)
colnames(samplematrix)<-indexer

ann.map.full<-do.dijkstra(map)

soln.path.full<-dijkstra.soln(ann.map.full)
dcost.full<-path.cost(soln.path.full)
allratios.full<-totalratio(soln.path.full)
finalratio.full<-prod(allratios.full)
freeen.full<-log(finalratio.full)

#############
## PART 3: Optimal path in restricted grid
#############

dof<-c(TRUE, FALSE)
source('varBAR_dijkstra_import.r')

ann.map.restrict<-do.dijkstra(map)

soln.path.restrict<-dijkstra.soln(ann.map.restrict)
dcost.restrict<-path.cost(soln.path.restrict)
allratios.restrict<-totalratio(soln.path.restrict)
finalratio.restrict<-prod(allratios.restrict)
freeen.restrict<-log(finalratio.restrict)

###############
## PART 4: Visualize paths - param view only
###############

full_df<-as.data.frame(do.call(rbind, soln.path.full))
full_df<-cbind(full_df, ind='Augmented')
restrict_df<-as.data.frame(do.call(rbind, soln.path.restrict))
restrict_df<-cbind(restrict_df, ind='Standard')
paramplot_df1<-rbind(full_df, restrict_df)
colnames(paramplot_df1)<-c('lambda','temp','ind')

p1<-ggplot(data=point.grid, aes(x=lambda, y=temp))
a1<-geom_path(data=paramplot_df1, aes(x=lambda, y=temp, colour=ind), arrow=arrow(), size=2)
a2<-theme(legend.position='bottom')
a3<-scale_colour_discrete(name='')

paramplot1<-p1+a1+a2+a3+coord_fixed(ratio=0.2) + geom_point()

ggsave(file='varBAR-uniuni-paramview.pdf', plot=paramplot1)
