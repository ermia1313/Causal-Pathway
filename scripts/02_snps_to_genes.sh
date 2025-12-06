# Running kirk.exe and spock.exe with work.txt

# Combining the snps_list_hugo.txt.* into list_hugo_CAUSE.txt:

echo "hugo wELPD wELPD^2 oELPD r2 tagged_SNP tag_SNP" \
  | tr ' ' $'\t' > ./list_hugo_CAUSE.txt

cat snps_list_hugo.txt.* \
  | sort -t $'\t' -k 2,2 -g \
  >> ./list_hugo_CAUSE.txt

# Any checks like grep -H -w DCC & CAMKV.

# Also move existing snps_to_genes.bat into this folder:

./snps_to_genes.bat


