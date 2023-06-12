# Zamborain-Masonetal_2023_ReefSustainability
This repository contains the code to implement the analyses from the manuscript titled "Sustainable reference points for multispecies coral reef fisheries" (Zamborain-Mason et al. 2023).
This project 
(i) estimates context-specific sustainable reference points multispecies reef fish assemblages based on environmental conditions 
(ii) assesses the status of multispecies reef (fish) fisheries using available estimates of catch and biomass
(iii) estimates the trade-offs between long-term fisheries production and other ecological fish metrics 

Methods are detailed in the manuscript and supplementary material.
If you use this data or code please cite the paper accordingly. 

# System requirements
All analyses were implemented in R (R version 4.2.1 (2022-06-23)). For the main analyses we used the RSTAN package. To account for methodological effects in our ecosystem response variables we used the brms package. Additional packages were used to re-arrange the data, perform sensitivity analyses, check model fit and visualization. Specific packages used can be found in code files (i.e.,  `ZamborainMason_ReefSustainability_REVISED.R`, `ZamborainMason_ReefSustainability_sp_specific_RESVISED.R`, and `ZamborainMason_ReefSustainability_workflow_REVISED.R`).

# Installation guide
If you haven't already, start by downloading R and RStudio on your computer. 

  - R is an open-source project of programming language, that is easily available for free download via the web. 
First, go to https://cran.r-project.org/ and click on the version of R that is suitable for your computer (e.g., Windows or Mac). 
If you pressed on the Mac link, press on the first link that appears after "Latest release:" (e.g., R-4.1.2.pkg). 
If you pressed on the Windows link, press the link names "base" just under "Subdirectories:". In the next page, press on the link to download the latest version fo R (e.g., Download R 4.1.2 for Windows). 
R will start downloading on to your computer.  Once it has finished downloading, open it and go through the installation process (.exe file in Windows). follow the default for your system type. For some, it will make you choose between 32 and 64-bit, which refers to the way a computer processor handles information.You can access your computer's system type by searching for "System information" in your windows search icon. 

  - RStudio is an integrated development environment that uses R (some classify it as a more "user friendly" way of using R). It is also easily available for free download via the web. 
To download RStudio, go to https://www.rstudio.com/products/rstudio/download/. Click the download button under the RStudio Desktop Open Source License (Free). Then, under "All installers", choose the link that matches your computer (e.g., Windows or Mac). 
RStudio will download to your computer. Open it and go through the installation process (.exe file for Windows). 

Install the required packages within R or RStudio and you are set-up to perform the analyses! For your code to run, make sure you have all data, code, and dependencies (e.g., see bellow) in the correct working directory

# Data and code
Main R code:
  - `ZamborainMason_ReefSustainability.R`: This commented code uses the main datasets and stan models to guide the user from top to bottom to perform the 
    main analyses and figures. It also includes major sensitivity analyses (e.g., choice of surplus production model or catch statistics).
    _Note that, given the nature of our models, results may vary slightly from run to run_
  
Main datasets:
  - `jurisdictionscale_data.csv`: cleaned jurisdiction-level data used in our main analyses.
  - `reefscale_data.csv`: cleaned reef site data used in our main analyses.

Stan models: 
  - `Null_gomperztfox_grav.stan`: Null reference point model.
  - `Full_gompertzfox_grav.stan`: Main reference point model (i.e., full).
  - `Full_gompertzfox_noremote_grav.stan`: reference point model not including remote locations (sensitivity analysis).
  - `Full_gompertzfox_exportrate.stan`: reference point model modelling exports as a rate (sensitivity analysis).
  - `Full_gompertzfox_dig.stan`: reference point model modelling exports as a proportion of the intrinsic growth rate (sensitivity analysis).
  - `Full_gompertzfox_dig_fixed.stan`: reference point model modelling exports fixed as a proportion of the intrinsic growth rate (sensitivity analysis).
  - `Full_pellatomlinson_free_grav.stan`: reference point model under a pella-tomlinson surplus production  (sensitivity analysis).
  - `Full_schaefer_grav.stan`: reference point model under a graham-schaefer surplus production  (sensitivity analysis).
  - `Full_pellatomlinson_fixed_grav.stan`: reference point model under a pella-tomlinson surplus production fixing parameter value (sensitivity analysis).

Additional datasets, dependencies or code: 
  - `individualspscale_data.csv`: species-abundance SERF data used in our supplementary analyses.
  -  `ZamborainMason_ReefSustainability_sp_specific.R`: code to perform species-specific sensitivity analyses.
  - `ZamborainMason_ReefSustainability_workflow.R`: code to perform bayesian workflow.
  - `monitornew_Vehtarietal2019.R`: functions from Vehtari et al 2019* to monitor convergence.Code can be accessed in https://github.com/avehtari/rhat_ess
  - `monitorplot_Vehtarietal2019.R`: functions from Vehtari et al 2019* to monitor convergence.Code can be accessed in (https://github.com/avehtari/rhat_ess)
     *Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, Paul-Christian B?rkner (2019): Rank-normalization, folding, and localization: An improved               Rhat for assessing convergence of MCMC. arXiv preprint arXiv:1903.08008.
  - `wormetal2009_californiacurrentsp.csv`: system, genus, and species information from the "California current" system extracted from the supplementary  
     information of  Worm et al. 2009 (https://www.science.org/doi/10.1126/science.1173146) and used for our supplementary analyses. If you use these please cite
     the paper accordingly. 
  - `wormetal2009_spdata.csv`: system, genus, and species information extracted from the supplementary  
     information of  Worm et al. 2009 (https://www.science.org/doi/10.1126/science.1173146) and used for our supplementary analyses. If you use these please cite
     the paper accordingly. 


