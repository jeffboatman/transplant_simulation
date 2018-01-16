Simulation code for estimating the anticipated survival for following the treatment
regime "decline all low-quality organs" assuming either 1) no other patients are following the 
regime, or 2) all patients are following the regime.

main_sim.R is the main simulation file with both estimators implemented.

sim_truth_all_follow.R is used to estimate the true survival distribution assuming
all patients are following the treatment regime.

sim_one_follows.R is used to estimate the true survival distribution assuming
no other patients are following the treatment regime.

summary.R reads the outfiles and creates the simulation table in the manuscript.

Each file except summary.R was written to be run on 5 student servers (carbon, cesium, chromium, potassium, silicon). 
These can be modified by the user as needed.
