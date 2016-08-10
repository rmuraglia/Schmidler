# sbar_v_ais_viz.r

options(warn=-1)
suppressMessages(library(ggplot2))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
options(warn=0)

dat<-read.table('sbar_v_ais_out.txt', header=F, sep='\t')

colnames(dat)<-c('rep', 'smc_meth', 'targ_dist', 'num_distr', 'path_height', 'numtraj', 'prc_F', 'nummetrostep', 'estimate')
levels(dat[,2])<-c('AIS', 'seqBAR')
protocol<-paste(dat[,2], '\n(', dat[,4], ', ', dat[,6], ', ', dat[,8], ')', sep='')
true_ratio<-11.20187 # to check, just load yacas_precompute and see Z1/Z0
fracerr<-(dat[,9]-true_ratio)/true_ratio
dat<-cbind(dat, fracerr, protocol) # columns 10 and 11

#########
# PANEL A: same comput, compare seqbar to ais: different outcomes
#########

subA<-subset(dat, dat[,4]==6 & dat[,6]==600 & dat[,8]==10)
subA[,11]<-droplevels(subA[,11])

baseA<-ggplot(data=subA, aes(x=protocol, y=fracerr)) + geom_boxplot() + geom_hline(yintercept=0, colour='red', linetype='dashed') + theme(axis.text.x=element_text(size=10, angle=0, hjust=0.5, colour='black'))
labelA1<-ylab(expression(frac(hat(r)-r, r)))
labelA2<-theme(axis.title.y=element_text(angle=0))
plotA<-baseA + labelA2 + labelA1 + coord_cartesian(ylim=c(-0.2, 0.15)) + labs(title='A') + theme(plot.title=element_text(hjust=0, size=20))

######
# PANEL B: compare computation required for same result
######

subB1<-subset(subA, subA[,2]=='seqBAR')
subB2<-subset(dat, dat[,2]=='AIS' & dat[,4]==11 & dat[,6]==600 & dat[,8]==10)
subB3<-subset(dat, dat[,2]=='AIS' & dat[,4]==6 & dat[,6]==3000 & dat[,8]==10)
subB4<-subset(dat, dat[,2]=='AIS' & dat[,4]==6 & dat[,6]==600 & dat[,8]==20)

subB<-rbind(subB1, subB2, subB3, subB4)
subB[,11]<-droplevels(subB[,11])

baseB<-ggplot(data=subB, aes(x=protocol, y=fracerr)) + geom_boxplot() + geom_hline(yintercept=0, colour='red', linetype='dashed') + theme(axis.text.x=element_text(size=10, angle=0, hjust=0.5, colour='black'))
labelB1<-ylab(expression(frac(hat(r)-r, r)))
labelB2<-theme(axis.title.y=element_text(angle=0))
plotB<-baseB + labelB2 + labelB1 + coord_cartesian(ylim=c(-0.2, 0.15)) + labs(title='B') + theme(plot.title=element_text(hjust=0, size=20))

outplot<-arrangeGrob(plotA, plotB, nrow=1, widths=c(2,3))

ggsave(file='seqBARvAIS-uni-bi-v02.png', width=8, height=4, dpi=600, plot=outplot)




