### Reproducibility Archive for paper "The Impact of Ordinal Scales on Gaussian Mixture Recovery"

The below describes how to use the code files in the reproducibility archive to fully reproduce all results in the paper.

#### Simulation Study

- `GenerateModels.R` generates the true GMMs used in the simulation study that vary across the number of components $K$, the pairwise KL-divergence $\text{D}_{\text{KL}}$ and the number of variables $p$
- `aux_functions.R` contains various helper functions for generating data, fitting GMMs, and plotting
- `Simulation_LISA.R` contains the simulation code to be executed on a single node of the UvA LISA cluster system. On each node, I parallelize across 12 cores using the `foreach` package. I run this script on $100$ nodes to obtain the $100$ simulation iterations reported in the paper. The output of each script is a `.RDS` file with the results of that iteration.
- The folders /output_EEE and /output_VVV each contain the $100$ output files from `Simulation_LISA.R` for the constrained and unconstrained estimation, respectively
- `Evaluation.R` takes the simulation results in one of the two folders as input (this can be selected at the beginning) and then preprocesses the results and plots all figures in the paper (and additional ones not shown in the paper)
- `submit_jobs.sh` and `submit_all.sh` are bash-scripts that are used to send `Simulation_LISA.R` to the nodes on the UvA LISA cluster system

The simulation on each node completed in less than 1h with 2.10 GHz cores.

#### Illustration Figures

- `Illustration_Figures.R` contains the code to produce Figures 1, 2 and 4 in the paper


#### Additional Analyses

- `KLD_vs_Difficulty.R` contains the code to produce the results in Appendix C


#### Session Info

This is the session info on the nodes used for the simulation:

> sessionInfo()
R version 4.0.2 (2020-06-22)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Debian GNU/Linux 10 (buster)

Matrix products: default
BLAS:   /sara/eb/AVX2/Debian10/EB_production/2020/software/R/4.0.2-intel-2020a/lib/R/lib/libR.so
LAPACK: /sara/eb/AVX2/Debian10/EB_production/2020/software/R/4.0.2-intel-2020a/lib/R/modules/lapack.so

locale:
 [1] LC_CTYPE=en_US       LC_NUMERIC=C         LC_TIME=en_US       
 [4] LC_COLLATE=en_US     LC_MONETARY=en_US    LC_MESSAGES=en_US   
 [7] LC_PAPER=en_US       LC_NAME=C            LC_ADDRESS=C        
[10] LC_TELEPHONE=C       LC_MEASUREMENT=en_US LC_IDENTIFICATION=C 

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
[1] mclust_5.4.6      MASS_7.3-51.6     doParallel_1.0.15 iterators_1.0.12 
[5] foreach_1.5.0    

loaded via a namespace (and not attached):
[1] compiler_4.0.2   codetools_0.2-16
