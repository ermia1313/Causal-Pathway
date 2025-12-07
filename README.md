# Causal Pathway

# Causal and Pleiotropic Biological Pathways

This repository contains code to identify **causal** and **pleiotropic** biological pathways by integrating CAUSE outputs with curated gene sets (GO, HPO, KEGG, Reactome, and others).  

Starting from SNP-level CAUSE results (ELPD deltas), we:
1. Propagate signal through LD to nearby SNPs and genes.
2. Aggregate gene-level evidence within pathways.
3. Use large-scale randomization to obtain empirical p-values.
4. Classify pathways as predominantly **causal** or **pleiotropic**.
5. Summarize and de-duplicate GO terms using semantic clustering.

---

## Conceptual Overview

### From CAUSE SNPs to gene-level scores

CAUSE assigns each tag SNP an **ELPD delta** value:

- **Negative delta** → consistent with **causation**  
- **Positive delta** → consistent with **pleiotropy**

Because CAUSE works on LD-independent **tag SNPs**, we use PLINK to recover their LD “buddies”:

- Reference: 1000 Genomes, CEU, hg19  
- PLINK options: `--show-tags`, `r² ≥ 0.1`, window 250 kb  

For each tag SNP:

1. We identify all **tagged SNPs** in its LD block (r² ≥ 0.1).
2. Each tagged SNP inherits the tag’s ELPD delta, **scaled by its LD**:
   \[
   \text{wELPD}_{\text{tagged}} = \text{ELPD}_\text{tag} \times r^2(\text{tagged}, \text{tag})
   \]
3. SNPs are intersected with gene coordinates; each gene inherits the wELPD values of intersecting SNPs.
4. For each gene, we retain the most extreme (largest |delta|) LD-weighted ELPD, so longer genes do not dominate purely by size (LD decays with distance).

This produces a file of the form:

- `hugo` – official gene symbol  
- `wELPD` – LD-weighted ELPD score  
- `wELPD²` – squared score  
- `oELPD` – original ELPD for the tag SNP  
- `r2` – LD r² between tagged SNP and tag SNP  
- `tagged_SNP` – SNP inside the LD block  
- `tag_SNP` – original tag SNP used by CAUSE  

---

## Pathway Enrichment With Empirical P-values

For each pathway (gene set), we test whether genes **inside** the pathway show more extreme ELPD deltas than genes **outside** the pathway.

We perform **two separate tests** per pathway:

- **Causation test** – using genes with **negative** wELPD
- **Pleiotropy test** – using genes with **positive** wELPD

Let:

- \( N \) = number of in-pathway genes considered (negative or positive, depending on the test)
- \( S_{\text{obs}} \) = sum of wELPD values for those \( N \) in-pathway genes
- \( S_{\text{rand}} \) = sum of wELPD values from \( N \) genes randomly drawn from the appropriate **outside-pathway** pool

We generate up to \( T = 10^7 \) random samples and compute an empirical p-value:

\[
p = \frac{1}{T} \sum_{i=1}^{T} I\left(\left|S_{\text{rand}, i}\right| \ge \left|S_{\text{obs}}\right|\right)
\]

Interpretation:

- If **\( S_{\text{rand}} > S_{\text{obs}} \)** is frequent → pathway behaves more like a **pleiotropic** set (positive direction).
- If **\( S_{\text{rand}} < S_{\text{obs}} \)** is frequent → pathway behaves more like a **causal** set (negative direction).

Multiple testing correction is performed within each direction:

- FDR adjustment using `p.adjust(..., method = "fdr")`
- Pathways are then ordered by p-value and effect size (|sum|).

---

## Gene Set Databases

We currently support multiple pathway collections (GMT files), including:

- DPs  
- DrugBank  
- GO Biological Process (GOBP)  
- GO Cellular Component (GOCC)  
- GO Molecular Function (GOMF)  
- Hallmark  
- HumanCyc  
- IOB  
- KEGG  
- MSigDB  
- NCI  
- Panther  
- PathBank  
- PFOCR  
- Reactome  
- TFs  
- WikiPathways  

These are provided as GMT files, e.g.:

```text
Human_DPs_fgsea.gmt
Human_DrugBank_fgsea.gmt
Human_GOBP_fgsea.gmt
Human_GOCC_fgsea.gmt
Human_GOMF_fgsea.gmt
Human_Hallmark_fgsea.gmt
Human_HumanCyc_fgsea.gmt
Human_IOB_fgsea.gmt
Human_KEGG_fgsea.gmt
Human_MSigDB_fgsea.gmt
Human_NCI_fgsea.gmt
Human_Panther_fgsea.gmt
Human_PathBank_fgsea.gmt
Human_PFOCR_fgsea.gmt
Human_Reactome_fgsea.gmt
Human_TFs_fgsea.gmt
Human_WikiPathways_fgsea.gmt
