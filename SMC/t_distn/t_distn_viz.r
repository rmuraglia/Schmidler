# t_distn_viz.r

library(ggplot2)

dat<-read.table('t_distn_out.txt', header=F, sep='\t')
colnames(dat)<-c('rep', 'smc_meth', 'targ_dist', 'num_distr', 'path_height', 'numtraj', 'prc_F', 'nummetrostep', 'estimate')

true_ratio<-1.253314 
fracerr<-(dat[,9]-true_ratio)/true_ratio
dat<-cbind(dat, fracerr)

p1<-ggplot(data=dat, aes(x=smc_meth, y=fracerr)) + geom_boxplot() + geom_hline(yintercept=0, colour='red', linetype='dashed')
p1 + coord_cartesian(ylim=c(-0.15, 0.15))

# ggsave(file='pcrooks.jpg')