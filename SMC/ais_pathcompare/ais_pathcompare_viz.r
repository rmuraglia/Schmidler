# ais_pathcompare_viz.r

options(warn=-1)
suppressMessages(library(ggplot2))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
options(warn=0)

dat<-read.table('ais_pathcompare_out.txt', header=F, sep='\t')

colnames(dat)<-c('rep', 'smc_meth', 'targ_dist', 'num_distr', 'path_height', 'numtraj', 'prc_F', 'nummetrostep', 'estimate', 'var')

combined_estim<-array(NA, dim=c(5, 2))
combined_estim[,1]<-c(0:4)
for (i in 1:5) { 
    subdat<-subset(dat, dat[,5]==(i-1))
    c_estim<-sum(subdat[,9]/subdat[,10])/sum(1/subdat[,10])
    combined_estim[i,2]<-c_estim
}
colnames(combined_estim)<-c('path_height', 'combined_estim')
combined_estim<-as.data.frame(combined_estim)

# right side panel: ais estims by path height
p1<-ggplot(data=dat, aes(x=path_height, y=estimate)) + geom_boxplot(outlier.size=0) + geom_point(position_jitter(w=0.1, h=0)) + geom_point(data=combined_estim, aes(y=combined_estim), colour='red', shape=18, size=4) + geom_hline(yintercept=1, colour='red', linetype='dashed', size=1.25) + labs(title='AIS ratio estimates')
p1_zoom<-p1 + coord_cartesian(ylim=c(0.65, 1.5))

# left isde panel: distribution path view
numdistr<-21
lambda_vals<-round(seq(from=0, to=1, length.out=numdistr), 3)

opt_T<-function(lambda, b) { round(sqrt(b^2 - b^2/0.5^2 * (lambda-0.5)^2) + 1, 3) }

for (i in 0:4) {
    tempdat<-cbind(i, lambda_vals, opt_T(lambda_vals, i))
    if (i==1) { distr_df<-as.data.frame(tempdat) } # if first one, make a dataframe for storage
    else { distr_df<-rbind(distr_df, tempdat) }  # otherwise just append to existing one
}
colnames(distr_df)<-c('path_height', 'lambda_vals', 'T_vals')
distr_df[,1]<-as.factor(distr_df[,1])

p4<-ggplot(data=distr_df, aes(x=lambda_vals, y=T_vals, colour=path_height)) + geom_point(size=3) + geom_line(size=1.5) + theme(legend.position='bottom') + labs(title='Overview of paths')

dev.new(width=10, height=5)
grid.arrange(p4, p1_zoom, nrow=1)
# dev.off()



