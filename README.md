# Plant-FATE runs for Amazon 

This file documents the model specification and model runs for the Amazon paper (Joshi et al, in prep) and Plant-FATE's contribution to Amazon FACE MIP (AmzMIP, paper led by Carolina Blanco).

## Model code

Calibrations (for the Plant-FATE paper) were performed on commit `4f9867e` (unpatched code). However, this commit did not have all the cohort properties, and did not write instantaneous output (all outputs were written annually).

Hence I patched `CommunityProperties` from commit `dmip_v1.0` to output inst values, and commented out soil related outputs. Later when merging `feature/soil` into develop, these soil-related changes will conflict, but the originals are preserved in block comments, which can be resolved during merge. This patched code is on the following commit, where the final simulations are performed.

```         
* 5ac79bd [hotfix/output_props] update cohort props header
```

This commit will be merged into `develop` and into `release/oligomorphic_paper` before submission and a tag will be created for the first Plant-FATE release from the merge.

## Plant-FATE runs for Amazon manuscript 

### Simulation setup

For the Plant-FATE paper, all variables have been calibrated from relevant data. The input data, config files, and outputs reside on `GECO-Workstation-2` in `/data/scratch/jaideep/AmzMIP_server`. The folder structure is as follows:

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

### Lifetime fitness simulations

Simulations for calculating fitness landscapes are done using the R version of the patched code. The R script that does the calculations and plots is in `Rscripts/plot_fitness_wdopt_from_dir.R`. There a slight discrepancy in the max fitness species and the evolutionarily emergent species. This is simply because of average z\* over 100 years is used for fitness calculations. There is no difference when z\* does not fluctuate, e.g., in the ELE_evol run from `pfate_output_evol_vs_comm` .

### Figures

Rscripts for generating figures are in `Rscripts/`. Each script generates one figure.

## Plant-FATE runs for Amazon MIP

The MIP runs differ slightly from the runs for the Plant-FATE manuscript.

- For MIP, the model is not calibrated specifically with site data. Instead literature-derived parameter values are used for almost all parameters. Others are guesstimates.

- MIP simulations are done at a monthly timestep (0.08333 years)

### Representation of diversity

Diversity is represented in terms of three traits: 

   - Wood density (WD)
   - Max height
   - Xylem P50

   The community-level trait distribution is currently assumed to be unimodal, but with a finite (interspecific) variance. Under high-diversity scenarios, this variance allows the community-level mean (CWM) trait values to change over time following fitness gradients. 
   
   In principle, this variance also introduces correction terms to the calculation of emergent ecosystem properties, but for the current simulations, we have kept the variance small enough for these corrections to be negligible. The variance only contributes to the evolution of traits by gradually ascending the fitness landscape.
   
   Furthermore, in principle, the curvature of the fitness lanscape also contributes to changes in variance over time. However, we have not implemented this in Plant-FATE yet, so trait variances are assumed to be constant. 

   Under low-diversity scenarios, traits are held constant at the mode of the aforementioned distribution (variance is set to 0, such that there is no trait evolution). In both scenarios, only the mode appears in the model outputs.
   
### Model run commands

The commands for running the model are of the form:

```         
nohup plantfate config_files/p_amz_mip_final_decalibrated.ini -19999.99 20000.01 > logs/p_hist-ele_evol_4.txt &
```

### Simulation setup

The input data, config files, and outputs reside on `GECO-Workstation-2` in `/data/scratch/jaideep/AmzMIP_server`. Outputs for all scenarios are in the folder `pfate_output_mip`. Folders are named as `AmzMIP_HIST_AMB/ELE_evol/ld_*ky_c2_rs0.04`.

Final data submitted for the MIP is formatted using the R script `AmzMIP_format_outputs.R`. 

