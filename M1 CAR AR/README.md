## Model 1: CAR-AR(1) Structure
This repository contains three R scripts for simulation studies (Scenarios 1–3) under different group testing protocols: Dorfman testing (DT) and array testing (AT) with unknown sensitivity and specificity, and individual testing (IT) with known sensitivity and specificity.

Each script includes both the data generation process and the Markov chain Monte Carlo (MCMC) procedure for the corresponding protocol:
1. ##M1_AT1.R## -- AT protocol,
2. ##M1_DT1.R## -- DT protocol,
3. ##M1_IT1.R## -- IT protocol. 

These scripts are designed for parallel computing in a SLURM cluster environment and require 28 CPU cores.
