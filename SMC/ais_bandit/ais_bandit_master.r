# ais_bandit_master.r

# collect data for AIS bandit figure (5.1)

# set fixed params

smc_meth<-'ais'
targ_dist<-'norm'
num_distr<-21
numtraj<-200
prc_F<-1
nummetrostep<-25
total_trial<-90

# get variable params
args<-commandArgs(trailingOnly=TRUE)
if (length(args)==2) {
    strat<-as.character(args[1])
    rep<-sprintf("%03d", as.numeric(args[2]))
} else { 
    print('improper number of command line arguments supplied')
    print('quitting R...')
    quit('no')
}

# import methods
source('../smc_samplers.r')
source('../yacas_precompute.r')

# set up trial tracker
path_info<-array(NA, dim=c(total_trial, 4))
colnames(path_info)<-c('path', 'estimate', 'var')

# helper function
weighted.var <- function(x, w, na.rm = FALSE) {
    if (na.rm) {
        w <- w[i <- !is.na(x)]
        x <- x[i]
    }
    sum.w <- sum(w)
    sum.w2 <- sum(w^2)
    mean.w <- sum(x * w) / sum(w)
    (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm = na.rm)
}

# for equal strat, just choose from this sequence
equal_paths<-c(rep(0, length.out=10), rep(2, length.out=10), rep(4, length.out=10), rep(0, length.out=20), rep(2, length.out=20), rep(4, length.out=20))

# collect baseline info for each path
for (i in 1:10) {

    # set path
    path_height<-0
    temp_points<-T_get(lambda_points, path_height)
    distr_mat<-rbind(lambda_points, temp_points)
    distr_path<-as.list(as.data.frame(distr_mat))
    distr_indexer<-sapply(distr_path, paste, collapse='_')

    # get estimate
    draws<-collect_noneq_traj(prc_F)
    estimate<-smc_estim(draws[[1]][[1]], draws[[2]][[1]], draws[[1]][[2]], draws[[2]][[2]])

    # store info
    path_info[i,]<-c(path_height, estimate[1], estimate[2])
}

for (i in 11:20) {

    # set path
    path_height<-2
    temp_points<-T_get(lambda_points, path_height)
    distr_mat<-rbind(lambda_points, temp_points)
    distr_path<-as.list(as.data.frame(distr_mat))
    distr_indexer<-sapply(distr_path, paste, collapse='_')

    # get estimate
    draws<-collect_noneq_traj(prc_F)
    estimate<-smc_estim(draws[[1]][[1]], draws[[2]][[1]], draws[[1]][[2]], draws[[2]][[2]])

    # store info
    path_info[i,]<-c(path_height, estimate[1], estimate[2])
}

for (i in 21:30) {

    # set path
    path_height<-4
    temp_points<-T_get(lambda_points, path_height)
    distr_mat<-rbind(lambda_points, temp_points)
    distr_path<-as.list(as.data.frame(distr_mat))
    distr_indexer<-sapply(distr_path, paste, collapse='_')

    # get estimate
    draws<-collect_noneq_traj(prc_F)
    estimate<-smc_estim(draws[[1]][[1]], draws[[2]][[1]], draws[[1]][[2]], draws[[2]][[2]])

    # store info
    path_info[i,]<-c(path_height, estimate[1], estimate[2])
}

for (i in 31:total_trial) {

    # for the remainder of the trials, set path height based on current info and strategy

    # get path height info
    path0<-subset(path_info, path_info[,1]==0)
    path2<-subset(path_info, path_info[,1]==2)
    path4<-subset(path_info, path_info[,1]==4)
    var0<-weighted.var(path0[,2], 1/path0[,3])
    var2<-weighted.var(path2[,2], 1/path2[,3])
    var4<-weighted.var(path4[,2], 1/path4[,3])
    
    # choose based on strategy
    if (strat=='greedy') { # choose minimum variance
        temp_mat<-matrix(c(0,2,4, var0, var2, var4), nrow=2, ncol=3, byrow=TRUE)
        minvar_path<-which.min(temp_mat[2,])
        path_height<-temp_mat[1,minvar_path]
    } else if (strat=='proptnl') { # choose inversely proportional to variance (low variance = chosen more often)
        temp1<-1/c(var0, var2, var4)
        temp2<-temp1/sum(temp1)
        path_height<-sample(c(0, 2, 4), size=1, prob=temp2)
    } else if (strat=='equal') { # choose such that every path is run an equal amount of times
        path_height<-equal_paths[i]
    } else {
        print('strat should be one of "greedy", "proptnl" or "equal"')
        print('quitting R...')
        quit('no')
    }

    temp_points<-T_get(lambda_points, path_height)
    distr_mat<-rbind(lambda_points, temp_points)
    distr_path<-as.list(as.data.frame(distr_mat))
    distr_indexer<-sapply(distr_path, paste, collapse='_')

    # get estimate
    draws<-collect_noneq_traj(prc_F)
    estimate<-smc_estim(draws[[1]][[1]], draws[[2]][[1]], draws[[1]][[2]], draws[[2]][[2]])

    # store info
    path_info[i,]<-c(path_height, estimate[1], estimate[2])
}

# calculate combined estimate and squared error
true_ratio<-1
combined_estim<-sum(path_info[,2]/path_info[,3]) / sum(1/path_info[,3])
sq_err<-(true_ratio-combined_estim)^2

# print out
bandit_out<-matrix(data=c(rep, strat, combined_estim, sq_err), nrow=1)
names(bandit_out)<-c('rep', 'strategy', 'estimate', 'sqerr')
write.table(bandit_out, file='ais_bandit_out.txt', append=TRUE, quote=F, sep='\t', col.names=F, row.names=F)


