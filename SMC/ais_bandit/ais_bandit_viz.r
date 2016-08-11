# ais_bandit_viz.r

options(warn=-1)
suppressMessages(library(ggplot2))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
options(warn=0)

dat<-read.table('ais_bandit_out.txt', header=F, sep='\t')

colnames(dat)<-c('rep', 'strategy', 'estimate', 'sqerr')

baseC<-ggplot(data=dat, aes(x=strategy, y=sqerr)) + geom_boxplot() + theme(axis.text.x=element_text(size=10, angle=0, hjust=0.5, colour='black'))
labelC1<-ylab(expression((hat(r)-r)^2))
labelC2<-theme(axis.title.y=element_text(angle=0))
# uncomment for left aligned panel title, like in phRMA figure
plotC<-baseC + labelC2 + labelC1 + coord_cartesian(ylim=c(0.02, 0.075)) #  + labs(title='C') + theme(plot.title=element_text(hjust=0, size=20))

# ggsave(file='ais-bandit.png', width=4, height=4, dpi=600, plot=plotC)