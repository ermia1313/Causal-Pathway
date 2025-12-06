This is where someone currently has that big loop with inline R that builds pway_2024jun_${pway}.txt.

#!/usr/bin/env bash
set -e

for pway in DPs DrugBank GOBP GOCC GOMF Hallmark HumanCyc IOB KEGG MSigDB NCI Panther PathBank PFOCR Reactome TFs WikiPathways; do
  nb_pways=$(cat data/reference/Human_${pway}_fgsea.gmt | wc -l)

  echo -e "pway\tonto\tdesc\tdir\tsum\tpval\tpadj\tl_in\tl_out\tl_pway\tidx\tleading_edge" \
    > data/results/pway_2024jun_${pway}.txt

  for ipway in $(seq 1 ${nb_pways}); do
    if [ -f data/results/pway_2024jun_${pway}_${ipway}.txt ]; then
      cat data/results/pway_2024jun_${pway}_${ipway}.txt
    fi
  done >> data/results/pway_2024jun_${pway}.txt

  Rscript R/combine_pway_tables.R ${pway}
done
05_build_ensg_hugo_snp_table.sh
Contains the big join command that builds list_ensg_hugo_snp_delta_delta2.txt:

bash
Copy code
LANG=en_EN join -t $'\t' -1 1 -2 2 \
  <( tail -n +2 data/reference/mart_Hsap_Hsap.tsv | cut -f1,3 | LANG=en_EN sort -k1,1 ) \
  <( LANG=en_EN join -t $'\t' -1 1 -2 1 \
       <( zcat data/reference/rsID_ENSG.txt.gz ) \
       <( tail -n +2 data/raw/list_rsID_CAUSE.txt | LANG=en_EN sort -k1,1 ) \
     | LANG=en_EN sort -k2,2 ) \
  | awk 'BEGIN {FS="\t";OFS="\t"} { if( length($2) > 1 ) print $0, sqrt(($4)^2) }' \
  | sort -t $'\t' -k 5,5 -g -r \
  | awk 'BEGIN {FS="\t"} ( !seen[ $2 ]++ )' \
  | uniq \
  > data/results/list_ensg_hugo_snp_delta_delta2.txt
