Note: this directory contains code that runs without issue, but it is unclear if the output is correct. This is the directory for an incomplete final test case.

Here the initial state is high temperature, zero external field (disordered spins), and the target state is low temperature (below the critial point) and zero external field (ordered spins with two metastable states - majority spin up or majority spin down). 
The reaction coordinates are the (inverse) temperature (beta) taking the role of lambda, and an external magnetic field taking the role of the auxilliary variable.
Normally, we would expect a phase transition with no magnetic field, creating a hard to cross barrier with high variance. With the addition of an external field, we were hoping to overcome this phase transition, by allowing random spins at high temperature to align to an external field before settling into one of the two state states.
This example remains incomplete, as we observed consistent errors in the ratio estimates, as they differed from what we expected based on Onsager's exact solution. Furthermore, enumeration of states and manual calculation of the partition function revealed its own discrepancy with Onsager's solution, making it unclear where the issues lied.
Due to lack of time, the troubleshooting was abandoned.
