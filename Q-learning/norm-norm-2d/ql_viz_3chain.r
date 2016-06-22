# ql_viz_3chain.r

library(ggplot2)
library(grid)
library(gridExtra)
library(gplots)

idnum<-sprintf('%02d', repnum)
if (temp_dof) {
    filename<-paste('normnorm-3chains-', idnum, sep='')
} else {
    filename<-paste('normnorm-3chains-restrict-', idnum, sep='')
}
save.image(paste(filename, '.RData', sep=''))

# info panel
emptydf<-data.frame(0,0)
colnames(emptydf)<-c('x','y')
labelplot<-ggplot(emptydf, aes(x=x, y=y)) + geom_blank() + theme_bw() + labs(x='', y='')
title_string<-paste('Q-learning norm-norm \n mu0=', mu0, ', sig0=', sig0, '\n mu1=', mu1, ', sig1=', sig1, '\n ', num_traj, ' particles/episode \n converged in 30+', nrow(tim[[4]]), ' episodes \n tolerance level=', composite_tol, sep='')
infopanel<-labelplot+geom_text(x=0, y=0, label=title_string, size=3) + theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +theme(axis.ticks = element_blank(), axis.text.y = element_blank())

# first chain
path1<-t(sapply(tim[[1]][[1]], split_func))
colnames(path1)<-c('lambda', 'temperature')
path1<-as.data.frame(path1)
p1<-ggplot(data=point_grid, aes(x=lambda, y=temperature)) + geom_point() + geom_path(data=path1, arrow=arrow(), size=3) + labs(title=paste('ratio=', round(tim[[2]][1],4), ', var=', round(tim[[3]][1],4), sep=''))

# second chain
path2<-t(sapply(tim[[1]][[2]], split_func))
colnames(path2)<-c('lambda', 'temperature')
path2<-as.data.frame(path2)
p2<-ggplot(data=point_grid, aes(x=lambda, y=temperature)) + geom_point() + geom_path(data=path2, arrow=arrow(), size=3) + labs(title=paste('ratio=', round(tim[[2]][2],4), ', var=', round(tim[[3]][2],4), sep=''))

# third chain
path3<-t(sapply(tim[[1]][[3]], split_func))
colnames(path3)<-c('lambda', 'temperature')
path3<-as.data.frame(path3)
p3<-ggplot(data=point_grid, aes(x=lambda, y=temperature)) + geom_point() + geom_path(data=path3, arrow=arrow(), size=3) + labs(title=paste('ratio=', round(tim[[2]][3],4), ', var=', round(tim[[3]][3],4), sep=''))

outplot<-arrangeGrob(infopanel, p1, p2, p3, nrow=2)

ggsave(file=paste(filename, '.png', sep=''), width=8, height=8, dpi=100, plot=outplot)

# png(filename=paste(filename, '.png', sep=''), width=800, height=800, res=100)
# outplot
# dev.off()