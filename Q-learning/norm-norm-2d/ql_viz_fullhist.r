# ql_viz_fullhist.r

library(ggplot2)
library(grid)
library(gridExtra)
library(gplots)

idnum<-sprintf('%02d', repnum)
filename<-paste('normnorm-monitor-', idnum, sep='')
# filename<-'normnorm-monitor-01'

# unique(tim[[1]])
unique_solns<-unique(tim[[1]])
num_solns<-length(unique_solns)

plot_df<-array(NA, dim=c(num_episode, 4))
colnames(plot_df)<-c('iter', 'ratio', 'var', 'soln')
plot_df[,1]<-c(1:num_episode)
plot_df[,2]<-tim[[2]]
plot_df[,3]<-tim[[3]]
for (i in 1:num_solns) {
    inds<-which(sapply(tim[[1]], identical, unique_solns[[i]]))
    plot_df[inds, 4]<-i
}
plot_df<-as.data.frame(plot_df)
plot_df[,4]<-as.factor(plot_df[,4])

# get first path with real solution
# plot_df
tmp_ind<-min(which(!is.na(plot_df[,2])))
first_path<-plot_df[tmp_ind, 4]

filler_cols<-rep('black', as.numeric(first_path)-1)
pcolors<-c(filler_cols, rich.colors(num_solns-length(filler_cols)))

# info panel
emptydf<-data.frame(0,0)
colnames(emptydf)<-c('x','y')
labelplot<-ggplot(emptydf, aes(x=x, y=y)) + geom_blank() + theme_bw() + labs(x='', y='')
title_string<-paste('Q-learning norm-norm \n mu0=', mu0, ', sig0=', sig0, '\n mu1=', mu1, ', sig1=', sig1, '\n ', num_episode, ' episodes \n', num_traj, ' particles/episode')
infopanel<-labelplot+geom_text(x=0, y=0, label=title_string, size=3) + theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +theme(axis.ticks = element_blank(), axis.text.y = element_blank())

# ratios overview
ratioplot<-ggplot(data=plot_df, aes(x=iter, y=ratio, colour=soln)) + geom_point() + theme(legend.position='none') + scale_colour_manual(values=pcolors, breaks=c(first_path:num_solns))

# vars overview
varplot<-ggplot(data=plot_df, aes(x=iter, y=var, colour=soln)) + geom_point() + theme(legend.position='none') + scale_colour_manual(values=pcolors, breaks=c(first_path:num_solns))

# path overviews - pick and choose
path_overviews<-vector('list', num_solns)

for (i in 1:num_solns) {
    path_overviews[[i]]<-ggplot(data=point_grid, aes(x=lambda, y=temperature)) + geom_point() + geom_path(data=as.data.frame(unique_solns[[i]]), arrow=arrow(), colour=pcolors[i], size=3) + labs(title=paste('soln ', i))
}

path_overviews[[num_solns+1]]<-infopanel
path_overviews[[num_solns+2]]<-ratioplot
path_overviews[[num_solns+3]]<-varplot
plot_inds<-c((num_solns+1):(num_solns+3), as.numeric(first_path):num_solns)

print(length(plot_inds))

plot_args<-c(path_overviews[plot_inds], list(ncol=4, widths=c(1,1,1,1)))

write.table(plot_df, file=paste(filename, '.txt', sep=''), quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)

png(filename=paste(filename, '.png', sep=''), width=1200, height=1200, res=100)
do.call(grid.arrange, plot_args)
dev.off()
