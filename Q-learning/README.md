Q-learning directory contains methods for carrying out approximate, inline and adaptive path searching to optimize paths for minimum variance ratio estimation.
In all cases, we use the nonequlibrium pCrooks estimator defined in the _SMC_ directory to obtain BAR ratio estimates and varBAR variance estimates to use as edge weights at low computational cost.

_norm-norm-2d_ directory contains the last pieces for the thesis. Search for paths between mean-shifted normal distributions in a 2d grid indexed by lambda and temperature. Master file can be interconverted to run either a fixed length Q-learning search (visualized with _ql\_viz\_fullhist.r_), or to run an adaptive version which monitors convergence and ends when a convergence criterion is met (visualized with _ql\_viz\_3chain.r_).

_ising-2d-ferro_ directory contains files for a second example (two-dimensional square Ising ferromagnet with periodic boundaries) which was not completed in time for the thesis and remains incomplete. 
