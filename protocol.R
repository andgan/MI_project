##################################################
# MYOCARDIAL INFARCTION ANALYSIS PROTOCOL
##################################################

##################################################
# PIPELINE DESCRIPTION 
##################################################

# association_analysis_pipeine.R -> pipleine for analysis of ICD10 codes both before and after baseline
# sensitivity_analysis_pipeline.R -> pipeline for analysis of ICD10 codes only before baseline

##################################################
# PIPELINE DESCRIPTION FOR ASSOCIATION ANALYSIS
##################################################

# PART 1: Load and format data
# Description: Load and format data correctly. 

# PART 2: Clean data based on exclusion criteria
# Description: Clean data based on various exclusion criteria. Exclude or remove certain ICD10 codes or individuals.
# Note: Change the dates for variables "study.start" and "study.end". 

# PART 3: Survival analysis for PRS
# Description: Calculate c-indices and hazard ratios with survival analysis for PRS. 

# PART 4: Survival analysis for ICD10 codes 
# Description: Calculate c-indices with survival analysis for ICD10 codes. 

# PART 5: Analysis of interaction for population subsets of ICD10 codes
# Description: Analyze the interaction between ICD10 codes and PRS in certain population subsets. 

# PART 6: Analysis of continuous PRS for population with and without previous relevant comorbidities
# Description: Calculate c-indices and hazard ratios for continuous PRS. 

# PART 7: Analysis of dichotomous PRS for population with and without previous relevant comorbidities
# Description: Calculate c-indices and hazard ratios for dichotomous PRS. 

# PART 8: ROC analysis of PRS with ICD10 codes and dichotomous PRS
# Description: Calculate AUCs, sensitivities, and specificites. 

##################################################
# RESULTS (16 files total)
##################################################

# association_prs.tsv 					(Part 3)
# association_analysis.tsv				(Part 4)
# association_subset_interaction.tsv	(Part 5)
# association_cont_prs.tsv				(Part 6)
# association_dich_prs.tsv				(Part 7)
# association_icd_sens_spec.tsv			(Part 8)
# association_prs_auc.tsv				(Part 8)
# association_prs_roc.tsv				(Part 8)

# sensitivity_prs.tsv 					(Part 3)
# sensitivity_analysis.tsv				(Part 4)
# sensitivity_subset_interaction.tsv	(Part 5)
# sensitivity_cont_prs.tsv				(Part 6)
# sensitivity_dich_prs.tsv				(Part 7)
# sensitivity_icd_sens_spec.tsv			(Part 8)
# sensitivity_prs_auc.tsv				(Part 8)
# sensitivity_prs_roc.tsv				(Part 8)