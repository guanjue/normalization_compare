##################################
script_folder='/Volumes/MAC_Data/data/labs/hardison_lab/normalization_compare/bin/'
##################################
	input_folder='input_data/'
	rsem_matrix='rsem_matrix/'
	output_figures='output_figures/'
	###### initiate folders
	### mkdir rsem_matrix module folder
	if [ -d "$rsem_matrix" ]; then  
    	rm -r $rsem_matrix
	fi
	mkdir $rsem_matrix

	### mkdir output_figures module folder
	if [ -d "$output_figures" ]; then  
    	rm -r $output_figures
	fi
	mkdir $output_figures

##################################
	### rsem_data_list.txt
	echo rsem_data_list.txt
	ls /Volumes/MAC_Data/data/labs/zhang_lab/vision_project/rna_seq_pred/rsem/*.rsem.*.bed > $input_folder'rsem_bed_list.txt'

	### get coding gene matrix
	time python $script_folder'rna_seq_matrix/merge_bed_matrix.py' -l $input_folder'anno2.txt' -i $input_folder'rsem_bed_list.txt' -t $input_folder'target_gene_class.txt' -u 500 -d 500 -o $rsem_matrix

	### quantile normalization
	time python $script_folder'rna_seq_matrix/quantile_normalization.py' -i $rsem_matrix'all_rna_sampe_target_class.txt' -o $rsem_matrix'all_rna_sampe_target_class_qn.txt'

	### plot hist before and after qn norm
	time Rscript $script_folder'figures/rna_od_qn_hist.R' $rsem_matrix'all_rna_sampe_target_class_qn.txt' $rsem_matrix'all_rna_sampe_target_class.txt' $output_figures'all_rna_sampe_target_class_qn_hist.pdf' $output_figures'all_rna_sampe_target_class_hist.pdf'
