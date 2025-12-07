# It contains the big join command that builds list_ensg_hugo_snp_delta_delta2.txt

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
