# Pathogenic Variants in all of us

Using the ACMG v3.3 list of pathogenic variants

After running files in src from 0 to 3. (2.5 and 3 won't run as executables for some reason but you can copy and paste code. I've no idea why they don't produce the same result)

This will produce a bunch of vcfs containing only pathogenic and likely pathogenic variants inside all the genes from the stated list across several batched vcfs. 

Then we combine them into a single VCF.


TODO:

## Where are the missing genes in our pathogenic variant list?
These are all the genes we do have:
```
ACTA2,APC,APOB,BMPR1A,C12orf43,CACNA1S,CCDC40,CEP85L,DSP,DSP-AS1,FAS,GAA,GLA,HFE,HFE-AS1,HNF1A,KCNH2,KREMEN1,LOC106560211,LOC126861339,LOC129994371,MYL2,NF2,OTC,PCSK9,PLN,PMS2,PRKAG2,PTEN,RET,RPE65,RPL36A-HNRNPH2,SDHD,TGFBR1,TNNI3,TNNT2,TRDN
```
There are 37. 26 overlap with the set of pathogenic variants (others are found because a variant might relate to both). 26 pathogenic genes we found are:
```
ACTA2,APC,APOB,BMPR1A,CACNA1S,DSP,GAA,GLA,HFE,HNF1A,KCNH2,MYL2,NF2,OTC,PCSK9,PLN,PMS2,PRKAG2,PTEN,RET,RPE65,SDHD,TGFBR1,TNNI3,TNNT2,TRDN
```
Leacing 100-26=74 genes we did not capture. Are the present at all. I'm going to check out original clinvar list before we overlapped with any data from AllOfUs to make sure we have variants where both these genes exist and there are pathogenic or likely pathogenic variables available. List of missing genes that I'm going to test with. There are some repeats so the following 68 genes are the ones we are missing.
```
ABCD1           ACTC1           ACVRL1          ATP7B           BAG3            BAG3            BRCA1           BRCA2           BTD              CALM1           CALM1           CALM2           CALM2           CALM3           CALM3           CASQ2           COL3A1           CYP27A1         DES             DES             DSC2            DSG2            ENG             FBN1            FLNC             FLNC            FLNC            KCNQ1           LDLR            LMNA            MAX             MEN1            MLH1             MSH2            MSH6            MUTYH           MYBPC3          MYH11           MYH7            MYH7            MYL3             PALB2           PKP2            RB1             RBM20           RYR1            RYR2            SCN5A           SCN5A            SCN5A           SDHAF2          SDHB            SDHC            SMAD3           SMAD4           SMAD4           STK11            TGFBR2          TMEM127         TMEM43          TNNC1           TP53            TPM1            TSC1            TSC2             TTN             TTR             VHL             WT1     
```
I double checked and these 2 lists together really are the complete list. Totaling 84 unique pathogenic genes. I used grep to segment so let's try the syntax on just the list of pathogenic variants.

Looking at just 1 example of missing variants. There are many pathogenic variants in the clinvar dataset associated with just ABCD1 for example.
```bash
bcftools view -H data/clinvar.vcf.gz | grep -F "ABCD1:" | grep -e "CLNSIG=Likely_pathogenic" -e "CLNSIG=Pathogenic" | wc -l
484
```

Now we should check whether these specific variants are appearing in the AllOfUs biobank vcf. This is more complicated because the AllOfUs biobank is sharded but we can find the ranges.
Range of ABCD1 pathogenic variants:
```
#CHROM  POS     ID
X       153725248       11317
X       153743577       2679143
```
This range is definitely in the interval_list files I used to segment so there aren't issues used in its generation.

I need to examine the interval lists in AllOfUs to find out which ones overlap
```
mapfile -t interval_files < <(gsutil -u "${GOOGLE_PROJECT}" ls "${CLINVAR_BASE_PATH}"*.interval_list 2>/dev/null)
INTERVAL_LISTS="output/interval_lists" && mkdir -p "$INTERVAL_LISTS"
# Stream input in because we have 20,000 interval_lists and can't add them all as arguments.
printf '%s\n' "${interval_files[@]}" | gsutil -u "${GOOGLE_PROJECT}" -m cp -I "$INTERVAL_LISTS/"
```
This copies all interval lists to a dir. Now I want to filter those interval lists for only those with overlap to `data/pathogenic_genes_1000000bp.interval_list` which we previously verified does have ABCD1 inside it:
```bash
cat data/pathogenic_genes_1000000bp.interval_list | grep ABCD1
chrX    152724495       154744755       +       ABCD1
```

Which interval files contain chrX
```bash
chrX_interval_lists=output/chrX/chrX_interval_lists
mkdir -p output/chrX
grep -l "^chrX" output/interval_lists/* > "${chrX_interval_lists}"
wc -l "${chrX_interval_lists}"
```
Then we must intersect. Can we feed in an intersect list of files?
```bash
GENE_INTERVAL_LIST="data/pathogenic_genes_1000000bp.interval_list"
overlap_file="output/chrX/overlap.interval_list"
chrX_input_tags="output/chrX/chrX_input_tags"
# Make input file
cat "${chrX_interval_lists}" | sed 's/^/-I /' > "${chrX_input_tags}"
gatk IntervalListTools \
        --arguments_file "${chrX_input_tags}" \
        -SI "${GENE_INTERVAL_LIST}" \
        --ACTION OVERLAPS \
        -O "${overlap_file}" \
        --QUIET false
```
The output interval_list file shows ranges overlapping ABCD1 min and max pathogenic variant. But we can't now go back to determine which of the interval_list files actually did this to figure out which vcf we should grab. 
Instead we need to iterate over all chrX interval_list files, overlap, if the overlap is not empty, then save that file name to an output file. Then we can use that file to grab the vcfs we are interested in.

We'll perform an overlap with every file and then just see if `grep "^chr"` returns anything because the header rows don't have this syntax.

```bash
intersected_interval_lists="output/chrX/intersected_interval_lists/"
intersected_interval_list_file="output/chrX/intersected_interval_paths.txt"
touch "${intersected_interval_list_file}"
mkdir -p "${intersected_interval_lists}"
mapfile -t lines < "${chrX_interval_lists}"
for line in "${lines[@]}"; do
  output_basename=$(basename "${line}")
  output_file="${intersected_interval_lists}""${output_basename}"".intersected"
  echo "${output_file}" >> "${intersected_interval_list_file}"
  gatk IntervalListTools \
    -I "$line" \
    -SI "${GENE_INTERVAL_LIST}" \
    --ACTION OVERLAPS \
    -O "${output_file}" \
    --QUIET true
done
```
Then we check for non-header output in these files. If it exists we know that that interval actual contains useful variants.
```bash
overlapped_interval_lists="output/chrX/overlapped_interval_lists.txt"
echo "" > "${overlapped_interval_lists}"
mapfile -t lines < "${intersected_interval_list_file}"
for line in "${lines[@]}"; do
  if grep -q -m 1 "^chr" "${line}"; then
    echo "${line}" >> "${overlapped_interval_lists}"
  fi
done
```
```bash
cat "${overlapped_interval_lists}"

0000019640.interval_list.intersected
to
0000019670.interval_list.intersected
-
0000019950.interval_list.intersected
to
0000019979.interval_list.intersected
-
0000020167.interval_list.intersected
to
0000020312.interval_list.intersected
```
This is 3 different ranges which is exactly what we wanted to see because we have 3 pathogenic genes on the x chromosome. And the last set of variants is exactly the lists that have ABCD1 which we lost. So this method will work.

So I repeated this whole process for all variants. 

Then I created parquet files from all vcf files to analyze in python using analysis/vcf_to_parquet.py

Then inside 



## Remaining Questions
- How many submitters do each pathogenic or likely pathogenic variant have in my set.
- What is the penetrance rate?
- What is the prevelance and penetrance by ancestry-related group
  - Do these relationships tend to follow patterns based on the type of variant?
  - Are there patterns in indvidual genes?
- What about variants in the region of a gene that may not be explicitly labeled with this gene
- Where are my missing genes?
  - The CALM gene family and a few others.


## Phenotyping
- 
- We have OMIM codes which we can convert to ORPHAcodes.
  - Some of these are 1-1 but most are not
- Then we have ORPHAcode 1-1 mappings to most SNOMED codes.
- This is all done using [ORPHA bindings](https://mappings.orphacode.org) manually online and in spreadsheets
Or
- Use [convert-pheno](https://github.com/CNAG-Biomedical-Informatics/convert-pheno) which can go between OMIM and OMOP
  - It's written in perl though. I've literally never used perl.
- [This paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10196319/) did OMIM > HPO > OMOP
  - Did they do it manually or are there nice tools?
Let's try convert-pheno 

## Collect ancestry information:
```
gsutil -u $GOOGLE_PROJECT cp gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv output/ancestry_preds.tsv
```


## Statistical Test:
- I need to more clearly define what I am comparing. 
- When is a statistical test valid if I very different case and control sizes?

Statistical test of 
- variants or 
- individuals 
Not sure what the difference is.

The more you can reproduce known results, then we can be more excited about novel results.
- If you generate a grouping. Are they really distinct?


The do associations with every variant with every phenotype.
Then:
- Is there a small set of variants that explains the total association? 
- Are there some variants that explain no association?

Then look for 2 instances of the same code for the same patient to call it a case.
Then we can look at pathogenecity

Phenotypic Analysis:
2. Parse presence of pathogenic variant in each individual not by clinvar label but by acmg definition of pathogenic which is slightly different for some genes. Then We can do summary statistics questions
3. Parse information about individuals
  - Phenotyping of relevant diseases
  - Ancestry
  - Covariates (age/sex/etc)

## What about familes?
Are certain prevalence? or prevalence w/o phenotype over represented within families?

Compare to known amounts otherwise noone will believe your novel results.


ACMG list from clingen vs from ACMG directly

gnomad vs clinvar variant list



exome sequencing may have higher coverage over rare variants.

plink file with everything.

survey data to compare individuals who have had a genetic test. 
- By race
- by depravation