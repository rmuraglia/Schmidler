# ising_import.r

# import methods for dealing with 2d ising ferromagnet sampling

################
# define energy functions
################

# overall unnormalized density function
# input: ONE config x as vector + temp and field params
# output: ONE scalar value for unnorm density
unnorm_dens<-function(x, beta, magn) {
    hamiltonian <- -coupl_energy(x, J) - magn_energy(x, magn)
    q_dens<-exp(-beta*hamiltonian)
    return(q_dens)
}

# try using a log density instead for more numerical stability
log_dens<-function(x, beta, magn) {
    hamiltonian <- -coupl_energy(x, J) - magn_energy(x, magn)
    log_out<- -beta*hamiltonian
    return(log_out)
}

# spin coupling contribution to energy
# for each pair of neighboring spins, calculate J*sig_i*sig_j
coupl_energy<-function(x, J) {
    total_en<-0
    for (i in 1:length(x)) {

        # get current spin value
        spin_i<-x[i]

        ## get neighbor spins 

        # if this spin is on the left edge, wrap around to get right edge spin
        if (i<=lattice_size) { 
            spin_L<-x[lattice_size*(lattice_size-1)+i]
        } else { 
            spin_L<-x[i-lattice_size]
        }

        # if spin is on the right edge, wrap to get left edge spin
        if (i>(lattice_size*(lattice_size-1))) {
            spin_R<-x[i-lattice_size*(lattice_size-1)]
        } else {
            spin_R<-x[i+lattice_size]
        }

        # if spin is on top
        if (i%%lattice_size==1) {
            spin_U<-x[i+lattice_size-1]
        } else {
            spin_U<-x[i-1]
        }

        # if spin is on bottom
        if (i%%lattice_size==0) {
            spin_D<-x[i-lattice_size+1]
        } else {
            spin_D<-x[i+1]
        }

        # get spin's contribution to coupling energy
        # divide by 2 to avoid double counting
        spin_en<-sum((J*spin_i/2) * c(spin_L, spin_R, spin_U, spin_D))
        total_en<-total_en+spin_en
    }
    return(total_en)
}

# external magnetic field contribution to energy
# for each spin, check if it is aligned: mu*sig_i
magn_energy<-function(x, mu) {
    spin_en<-x*mu
    total_mag_en<-sum(spin_en)
    return(total_mag_en)
}

#############
# sampling methods
#############

# note matrix command fills by column unless byrow argument used
# 11, 21, 31, 12, 22, 32, 13, 23, 33
# input: number of configurations to create (how many full lattices) + number of states along one edge of the square lattice (lattice is nxn)
# output: num_traj x n^2 matrix of spin values. one row = 1 lattice
generate_ising_hightemp<-function(num_traj, lattice_size) {
    # at high temp, just draw n^2 spins from c(-1, 1) randomly to create one starting configuration
    # do this num_traj times
    draws<-replicate(num_traj, sample(c(-1, 1), size=lattice_size^2, replace=TRUE))
    return(t(draws)) #transpose to get dimensions we want (1 row = 1 config)
}

generate_ising_lowtemp<-function(num_traj, lattice_size, majority_spin) {
    # at low temp, soln from paper indicates mean magnetization. unlikely for there to be even one spin that goes against the grain, but allow it anyway.
    if (majority_spin==1) {
        draws<-replicate(num_traj, sample(c(-1, 1), size=lattice_size^2, replace=TRUE, prob=c(0.000725, 0.999275)))
    } else if (majority_spin==-1) {
        draws<-replicate(num_traj, sample(c(1, -1), size=lattice_size^2, replace=TRUE, prob=c(0.000725, 0.999275)))
    } else {
        print('WARNING: for low temp, should have majority_spin flag set to +/- 1. \n Will generate random spins for this call.')
        draws<-t(generate_ising_hightemp(num_traj, lattice_size))
    }
    return(t(draws))
}

# note that random spin flip proposal density is symmetric, so we don't need that correction
ising_transition<-function(x, indexer_next) {
    prev_config<-x 
    next_state<-as.numeric(unlist(strsplit(indexer_next, split='_')))
    prev_log<-log_dens(prev_config, next_state[1], next_state[2])

    for (i in 1:num_metro_step) { # carry out a given number of metropolis steps per transition
        # select which spin(s) to flip
        flipped<-sample(c(1:lattice_size^2), size=num_spin_flip)

        # change those spin values
        trial_config<-prev_config
        trial_config[flipped]<-trial_config[flipped]*(-1)

        # get current density and trial density
        # prev_dens<-unnorm_dens(prev_config, next_state[1], next_state[2])
        # trial_dens<-unnorm_dens(trial_config, next_state[1], next_state[2])

        # log density of trial draw
        trial_log<-log_dens(trial_config, next_state[1], next_state[2])

        # get acceptance probability
        accept_prob<-exp(trial_log-prev_log)

        # accept/reject step, update draw?
        prob_draw<-runif(n=1, min=0, max=1)
        if (prob_draw<=accept_prob) {
            prev_config<-trial_config
            prev_log<-trial_log
        } 
    }
    return(prev_config)
}




