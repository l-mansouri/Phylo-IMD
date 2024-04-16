sed -i 's/exp_pdb_results_SUBSET49/mTMalign/g' reference_branches_titration_ME+3d+ML.R


for fam in `cat source_data/list_of_families_for_titration`; do
 Rscript reference_branches_titration_ME+3d+ML.R $fam ; 
 echo $fam
done

sed -i 's/mTMalign/exp_pdb_results_SUBSET49/g' reference_branches_titration_ME+3d+ML.R


python compute_titration_reference_branches.py source_data/list_of_families_for_titration exp
