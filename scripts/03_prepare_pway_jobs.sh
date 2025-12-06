# The loop

for pway in DPs DrugBank GOBP GOCC GOMF Hallmark HumanCyc IOB KEGG MSigDB NCI Panther PathBank PFOCR Reactome TFs WikiPathways; do
  nb_pways=`cat ./data/reference/Human_${pway}_fgsea.gmt | wc -l`
  for ipway in `seq 1 ${nb_pways}`; do
    echo /home/.../pway_empirical_P.bat ${pway} ${ipway}
  done
done | shuf > ./work.txt

# The kirk.exe + spock.exe calls for the pathway jobs.

# Move pway_empirical_P.bat into scripts/ too
