# ais_pathcompare_master.r

# collect data for AIS path comparison (fig 4.2)


# grab command line args
args<-commandArgs(trailingOnly=TRUE)

if (length(args)==8) {
    smc_meth<-as.character(args[1])
    targ_dist<-as.character(args[2])
    num_distr<-as.numeric(args[3])
    path_height<-as.numeric(args[4])
    numtraj<-as.numeric(args[5])
    prc_F<-as.numeric(args[6])
    nummetrostep<-as.numeric(args[7])
    rep<-sprintf("%03d", as.numeric(args[8]))
} else {
    print('improper number of command line arguments supplied')
    print('quitting R...')
    quit('no')
}

# import methods
source('../smc_samplers.r')
source('../yacas_precompute.r')

# collect draws
draws<-collect_noneq_traj(prc_F)

# get estimate
estimate<-smc_estim(draws[[1]][[1]], draws[[2]][[1]], draws[[1]][[2]], draws[[2]][[2]])

# set up print out info
print_out<-matrix(data=c(rep, smc_meth, targ_dist, num_distr, path_height, numtraj, prc_F, nummetrostep, estimate[1], estimate[2]))

# write to file
write.table(print_out, file='ais_pathcompare_out.txt', append=TRUE, quote=F, sep='\t', col.names=F, row.names=F)

