## Nemo + STRUCTURE pipeline
1. Simulation: `bash simulation/parallel_job_*` to run `./nemo` on a list of `ini` files
2. Fstat formatting: `bash fstat_format_rep_new.sh`
3. Fstat2str `fstat2str_rep.R`
4. Structure `bash structure/structure_threader_cor*.sh`
5. Merge Q across K `merge_Q_acorss_K_str_repeats.R` - generating `structure_rep_sum/*_merged_Q.txt`
6. Tidy Q across simulation `tidy_Q_across_nemo_simulation_rep.R` - generating
   1. `nemo_rep_sum/*_tidy_Q.tsv`: AQ_assignment_rate for each individual in each run of simulation - for individual STRUCTURE plot patterns and Q value distribution
   2. `nemo_rep_sum/*_sum_by_rep_pop_Q.tsv`: mean, var, freq_unadmixed for each simulation run (averaged across individuals within each population) - for summary Q values
7. Summarize and visualize results `summary_thresholdquantile90.R`

