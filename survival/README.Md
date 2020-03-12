#### Suvival Analysis results


Survival analysis carried out in the METABRIC cohort - breast cancer patients with normalized gene expression + metadata

Signatures obtained from HCC_MK and MDA data-sets are scored for each patient by summing the normalized expression value for genes in a given signature. 
Patients are then sorted into two groups - high score or low score - by setting a threshold. The threshold is identified as the value which maximizes the difference between two groups.
This can sometimes result in unbalanced groups (i.e 5 in group, 1995 in the other). 

For genes that are 'downregulated' --> "DOWN", expression is summed without multiplying by -1. 

Kaplan-Meir survival curves (km_plots) as well as hazard ratio plots (hazard_plots) are provided for each signature. 

The CSV file survival_results.csv consists of the following:
- pval_KM : The Kaplan-Meir p value - if signficant indicates signficant difference in # of people surviving in each group
- HZ_KM : the hazard ratio for being in one group vs othe other -- the larger the more signficant the effect, in either direction (> 1 or < 1)
--- Hazard ratio for various features such as whether the patient received Chemotherapy, and their PAM50 subtype
