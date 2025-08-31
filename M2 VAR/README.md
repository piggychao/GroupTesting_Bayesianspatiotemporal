## Model 2: VAR Structure  

This repository contains three R scripts for simulation studies (Scenarios 1–3) under different group testing protocols: Dorfman testing (DT) and array testing (AT) with unknown sensitivity and specificity, and individual testing (IT) with known sensitivity and specificity.  

Each script includes both the data generation process and the Markov chain Monte Carlo (MCMC) procedure for the corresponding protocol:  

1. **M2_AT1.R** – Array testing (AT) protocol  
2. **M2_DT1.R** – Dorfman testing (DT) protocol  
3. **M2_IT1.R** – Individual testing (IT) protocol  

These scripts are designed for parallel computing in a SLURM cluster environment and require 28 CPU cores.  
