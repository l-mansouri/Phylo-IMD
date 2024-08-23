for family_list in `ls source_data/titration_deciles`; do
 echo $family_list
 python compute_titration_reference_branches_quantiles.py source_data/titration_deciles/$family_list /home/luisasantus/Desktop/crg_cluster/newphylo/NF_TMalign/mTMalign/titration_every_5
done

