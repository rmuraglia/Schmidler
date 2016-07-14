# ql_master.r

# master file controlling Q-learning search for 2D Ising ferromagnet

# initial distribution: (0, 0) 
# infinite initial temperature (inverse temp is 0), and no external field
# spins are disordered and decoupled due to high temp. Configurations are uniformally distributed
# can use this as a reference partition function. For infinitely high temperature, Z=2^n. 2 spin choices per site with n sites.

# target distribtion: (1, 0)
# target beta is chosen to be 1 to match the PRE paper.
# with coupling constant J=1, we find that the critical temperature is b_crit = ln(1+sqrt(2))/2 ~= 0.4407. 
# choose beta_final = 1 ensures that we must cross this critical temperature barrier.
# at low temp, spins are tightly coupled. essentially only two valid configs: all spin up or all spin down. equilibrium per spin magnetization is m_eq = +/- 0.999275
# want to find partition function of this state.

# episodes are done as follows:
# initialize high temp state (draw from uniform distn)
# propagate towards low temperature state by sequential moves
# once reach target, initialize low temp state (draw each spin from 0.999275 vs 0.000725 weighted distn)
# propagate back towards high temp state by sequential moves
# calculate varBAR on this particle set.
# update R map, update Q map

print('Begin R script')
print(date())

# set working directory
switch(Sys.info()[['sysname']],
    Windows = {setwd('D:/work-duke/Schmidler/Q-learning/ising-2d-ferro')},
    Darwin = {setwd('/Users/rmuraglia/GitHub/Schmidler/Q-learning/ising-2d-ferro/')},
    Linux = {setwd('/dscrhome/rm204/ising-ql/')}
)

print('Loading required libraries...')
library(abind)

print ('Reading parameters from command line if present...')
args<-commandArgs(trailingOnly=TRUE)
if (length(args)==8) { # if appropriate number, use
    # repnum<-as.numeric(args[1])
    # magn_dof<-as.logical(args[2])
    # composite_tol<-as.numeric(args[3])
    move_jump<-as.numeric(args[1])
    beta_numpoints<-as.numeric(args[2])
    magn_numpoints<-as.numeric(args[3])
    magn_dof<-as.logical(args[4])
    num_traj<-as.numeric(args[5])
    num_episode<-as.numeric(args[6])
    num_metro_step<-as.numeric(args[7])
    rep<-as.numeric(args[8])
} else { # else use following defaults
    move_jump<-1
    beta_numpoints<-6
    # magn_numpoints<-9
    magn_numpoints<-5
    magn_dof<-FALSE
    num_traj<-100
    num_episode<-300
    num_metro_step<-25
    rep<-50
}

# move_jump<-1
# beta_numpoints<-6
# magn_numpoints<-9
# num_traj<-100
# num_episode<-600

# magn_dof<-TRUE
# magn_dof<-FALSE
# repnum<-01

# Ising model params
J<-1
beta_min<-0
beta_max<-1
magn_min<-0
magn_max_abs<-10
lattice_size<-10

# params for sequential sampling 
# num_metro_step<-25
num_spin_flip<-1

# params for all types of QL search
alpha<-0.8 # learning rate
gamma<-1 # discount factor

# params for full history QL search (fixed duration)
epsilon_init<-0
epsilon_tau<-floor(0.85*num_episode) # determines episilon rate of change - end searching on epsilon_tau'th search episode. after that, always choose lowest Q value move.

# params for convergence QL search
min_episode<-30
max_episode<-1000
# composite_tol<-0.025
# composite_tol<-0.05
chain_sd_mult<-1

source('ql_grid_setup.r')
source('pcrooks_import.r')
source('ising_import.r')
source('ql_algorithm.r')

# epsilon_tau<-2
# num_episode<-10
ql_out<-ql_full_history(q_map, r_map, alpha, gamma, epsilon_init, epsilon_tau)

if (magn_dof) {
    save.image(paste('ising-fullhist-psearch-TEST', sprintf('%03d', rep), '.RData', sep=''))
} else {
    save.image(paste('ising-fullhist-restrict-TEST', sprintf('%03d', rep), '.RData', sep=''))
}


# tim<-ql_episode(q_map, r_map, epsilon_init)
