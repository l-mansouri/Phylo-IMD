for fam in `cat source_data/list_families_fully_retrieved_uniprot_for_plot`; do
 Rscript reference_branches_titration_ME+3d+ML.R $fam ; 
 echo $fam
done

sed -i 's/exp_pdb_results_SUBSET49/af2_pred_pdb_results/g' reference_branches_titration_ME+3d+ML.R

for fam in `cat source_data/list_families_fully_retrieved_uniprot_for_plot`; do
 Rscript reference_branches_titration_ME+3d+ML.R $fam ; 
 echo $fam
done

sed -i 's/af2_pred_pdb_results/exp_pdb_results_SUBSET49/g' reference_branches_titration_ME+3d+ML.R


python compute_titration_reference_branches.py source_data/list_families_fully_retrieved_uniprot_for_plot AF2
python compute_titration_reference_branches.py source_data/list_families_fully_retrieved_uniprot_for_plot exp
