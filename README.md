# Plant-FATE runs for Amazon

This file documents the model specification and model runs for the Amazon paper (Joshi et al, in prep) and Plant-FATE's contribution to Amazon FACE MIP (AmzMIP, paper led by Carolina Blanco).

## Model code

Calibrations were performed on commit `4f9867e` (unpatched code). However, this commit did not have all the cohort properties, and did not write instantaneous output (all outputs were written annually).

Hence I patched `CommunityProperties` from commit `dmip_v1.0` to output inst values, and commented out soil related outputs. Later when merging `feature/soil` into develop, these soil-related changes will conflict, but the originals are preserved in block comments, which can be resolved during merge. This patched code is on the following commit, where the final simulations are performed.

```         
* 5ac79bd [hotfix/output_props] update cohort props header
```

This commit will be merged into `develop` and into `release/oligomorphic_paper` before submission and a tag will be created for the first Plant-FATE release from the merge.

## Simulation setup

For Amazon paper, all variables have been calibrated from relevant data. The input data, config files, and outputs reside on `GECO-Workstation-2` in `/data/scratch/jaideep/AmzMIP_server`. The folder structure is as follows:

-   `pfate_output_calib`: Contains outputs from various calibration trials, mainly calibrating `rs`. The final version is in `AmzMIP_HIST_AMB_evol_20ky_4`. Results are from the unpatched code.
-   `pfate_output_evol_vs_comm`: Contains 2 community simulations + 1 evolutionary simulation. Community runs are performed with 50 species for 20000 years each. Species fall along a gradient of wood density, and are identical in all other traits. One is run under aCO2 (414 ppm) and the other under eCO2 (828 ppm). The evolutionary run has a single species with evolvable wood density, and runs from -20k to 2020 under aCO2 and from 2020 to 20k under eCO2. These runs are used to validate the predictions of the evolutionary algorithm. These simulations are also performed with the unpatched code.
-   `pfate_output_co2scan`: Contains a scan of 20 runs with different eCO2 values, starting 414.2 ppm with increments of 50 ppm. Each run follows the MIP protocol until 2020 (368.9 ppm until 2000, historical CO2 from 2001-2020, and eCO2 from 2021 onwards). These form the bulk of results used for figures in Joshi et al. These runs use the patched code.
-   `pfate_output_supplementary`: Contains 2 simulations starting at 2000 under eCO2 = 614.2 ppm. One has +50% rs, the second has +50% zeta.

### Model run commands

The commands for running the model are of the form:

```         
nohup plantfate config_files/p_amz_mip_final.ini -19999.99 20000.01 > logs/p_hist-414_evol_4.txt &
```

They can be found in `run_commands.txt` .

## Lifetime fitness simulations

Simulations for calculating fitness landscapes are done using the R version of the patched code. The R script that does the calculations and plots is in `Rscripts/plot_fitness_wdopt_from_dir.R`. There a slight discrepancy in the max fitness species and the evolutionarily emergent species. This is simply because of average z\* over 100 years is used for fitness calculations. There is no difference when z\* does not fluctuate, e.g., in the ELE_evol run from `pfate_output_evol_vs_comm` .

## Figures

Rscripts for generating figures are in `Rscripts/`. Each script generates one figure.
