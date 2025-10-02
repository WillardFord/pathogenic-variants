# Pathogenic Variants in all of us

Using the ACMG v3.3 list of pathogenic variants

After running files in src from 0 to 3. (2.5 and 3 won't run as executables for some reason but you can copy and paste code. I've no idea why they don't produce the same result)

This will produce a bunch of vcfs containing only pathogenic and likely pathogenic variants inside all the genes from the stated list across several batched vcfs. 

Next steps:
1. Combine them into a single vcf
2. Parse presence of pathogenic variant in each individual not by clinvar label but by acmg definition of pathogenic which is slightly different for some genes. Then We can do summary statistics questions
3. Parse information about individuals
  - Phenotyping of relevant diseases
  - Ancestry
  - Covariates (age/sex/etc)