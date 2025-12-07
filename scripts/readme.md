# Pipeline Overview (Scripted Workflow)

The core pipeline is implemented as shell and R scripts. At a high level:

1. Input preparation:

* Start from list_rsID_CAUSE.txt (SNP + Delta_ELPD from CAUSE).

Use PLINK (--show-tags) to identify LD buddies for each CAUSE tag SNP (r² ≥ 0.1, 250 kb window, CEU reference, hg19).

2. SNP → gene mapping

Split a large SNP–gene mapping file (rsID_ENSG.txt.gz) into smaller files in shared memory.

For each block, call snps_to_genes.bat via a job list (work.txt) processed by kirk.exe and multiple parallel spock.exe workers.

Combine the per-block outputs into list_hugo_CAUSE.txt, containing genes and their LD-weighted ELPD scores.

3. Empirical pathway tests

For each pathway database (DPs, DrugBank, GOBP, GOCC, GOMF, Hallmark, HumanCyc, IOB, KEGG, MSigDB, NCI, Panther, PathBank, PFOCR, Reactome, TFs, WikiPathways):

Build a work list of jobs pway_empirical_P.bat <pway> <pathway_index>.

Run these jobs in parallel with kirk.exe + multiple spock.exe workers.

Combine per-pathway results into a single file pway_2024jun_<pway>.txt.

Apply FDR correction separately for causal (negative) and pleiotropic (positive) directions, and sort by p-value and effect size.

4. Summary and QC

Inspect top pathways per collection and direction (cause vs pleio).

Build list_ensg_hugo_snp_delta_delta2.txt, which links ENSG IDs, HUGO symbols, SNPs, and CAUSE deltas, retaining the most extreme delta per gene.

Optionally generate ECDF plots and other diagnostics for wELPD distributions.

Inputs and Outputs

# Key inputs

list_rsID_CAUSE.txt
SNPs and CAUSE ELPD deltas.

# SNP–gene mappings and annotations:

rsID_ENSG.txt.gz

mart_Hsap_Hsap.tsv (Ensembl ↔ HUGO mapping)

# Pathway / gene set GMTs:

Human_<DB>_fgsea.gmt (see list above)

# Key outputs

list_hugo_CAUSE.txt
Gene-level LD-weighted ELPD scores.
pway_2024jun_<DB>.txt
Empirical causal/pleiotropic pathway statistics and adjusted p-values for each database.
list_ensg_hugo_snp_delta_delta2.txt

Integrated table ENSG ↔ HUGO ↔ SNP ↔ CAUSE delta (most extreme per gene).

# Dependencies

1. PLINK (v1.9 or later) for LD calculations
2. R (with at least the following packages):
3. fgsea
4. foreach
5. doMC
6. rrvgo
7. A Unix-like environment with:
8. bash, awk, sort, join, zcat, etc.
9. Local compiled binaries:
kirk.exe, spock.exe (used as a simple parallel job dispatcher)
