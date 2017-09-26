import os
import numpy as np
################################################################################################
### read 2d array
def read2d_array(filename,dtype_used, sep):
	import numpy as np
	data=open(filename,'r')
	data0=[]
	for records in data:
		tmp = [x.strip() for x in records.split(sep)]
		data0.append(tmp)
	data0 = np.array(data0,dtype=dtype_used)
	data.close()
	return data0

################################################################################################
### write 2d matrix
def write2d_array(array,output):
	r1=open(output,'w')
	for records in array:
		for i in range(0,len(records)-1):
			r1.write(str(records[i])+'\t')
		r1.write(str(records[len(records)-1])+'\n')
	r1.close()

################################################################################################
def merge_bed_matrix(interval_and_label, rna_list, target_label, exp_up, exp_down, outputfolder):
	###### read inputs
	### read interval_and_label
	interval_and_label = read2d_array(interval_and_label, str, ' ')
	### read target_label
	target_label = read2d_array(target_label, str, '\t')[0,:]
	### read rna_list
	rna_list = read2d_array(rna_list, str, '\t')

	### read each rna sample
	all_rna_sampe = {}
	all_col_name = []
	for rna_sampel in rna_list:
		### get sample file location/name
		file_name = rna_sampel[0]
		### get sample file name => colnames
		col_name = file_name.split('/')[-1]

		### read sample file
		rna_sample = read2d_array(file_name, str, '\t')[:,:]

		### put sample file into dict
		all_rna_sampe[col_name] = rna_sample
		### get sample file order
		all_col_name.append(col_name)

	######
	### get target gene class positions 3rd col
	gene_class = interval_and_label[:,3]

	### initialize position vector
	indices = np.repeat(0, gene_class.shape[0])

	### include all target gene classes
	for c in range(0, target_label.shape[0]):
		indices = indices + (gene_class == target_label[c])

	### extract all target gene classes row id
	keep_position = np.where(indices != 0)

	### get all intervals
	interval_and_label = interval_and_label[keep_position, :][0,:,:]
	print('interval_and_label.shape')	
	print(interval_and_label.shape)
	### get strand info
	print(all_rna_sampe[all_col_name[0]].shape)
	strand = all_rna_sampe[all_col_name[0]][keep_position, 5]
	strand = np.transpose(strand)
	print('strand.shape')	
	print(strand.shape)

	### get bed files
	bed_matrix = []
	tss_bed_matrix = []
	for info, s in zip(interval_and_label, strand):
		chrom = info[0]	
		start = info[1]
		end = info[1]

		### get tss	
		if s[0] == '+':
			TSS = info[1]
			TSS_upstream = int(TSS) - exp_up
			TSS_downstream = int(TSS) + exp_down
		else:
			TSS = info[2]
			TSS_upstream = int(TSS) - exp_down
			TSS_downstream = int(TSS) + exp_up
		### get each bed vector
		bed_vector = [chrom, start, end, '0', '0', s[0]]
		tss_bed_vector = [chrom, TSS_upstream, TSS_downstream, '0', '0', s[0]]

		###
		bed_matrix.append(bed_vector)
		tss_bed_matrix.append(tss_bed_vector)



	all_rna_sampe_target_class = []
	for sample in all_col_name:
		all_rna_sampe_target_class.append(all_rna_sampe[sample][keep_position, 3])

	all_rna_sampe_target_class = np.transpose( np.array(all_rna_sampe_target_class)[:,0,:] )
	print(all_rna_sampe_target_class.shape)

	write2d_array(bed_matrix, outputfolder+'bed_matrix.bed')
	write2d_array(tss_bed_matrix, outputfolder+'tss_bed_matrix.bed')
	write2d_array(all_rna_sampe_target_class, outputfolder+'all_rna_sampe_target_class.txt')


############################################################################
############################################################################
#time python merge_bed_matrix.py -t anno2.txt -i rsem_bed_list.txt -t target_gene_class.txt -u 500 -d 500 -o outputfolder

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hl:i:t:u:d:o:")
	except getopt.GetoptError:
		print 'time python merge_bed_matrix.py -t interval_and_label -i rna_list -t target_label -u upstream_expand -d downstream_expand -o outputfolder'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python merge_bed_matrix.py -t interval_and_label -i rna_list -t target_label -u upstream_expand -d downstream_expand -o outputfolder'
			sys.exit()
		elif opt=="-l":
			interval_and_label=str(arg.strip())
		elif opt=="-i":
			rna_list=str(arg.strip())
		elif opt=="-t":
			target_label=str(arg.strip())		
		elif opt=="-u":
			exp_up=int(arg.strip())	
		elif opt=="-d":
			exp_down=int(arg.strip())	
		elif opt=="-o":
			outputfolder=str(arg.strip())		

	merge_bed_matrix(interval_and_label, rna_list, target_label, exp_up, exp_down, outputfolder)

if __name__=="__main__":
	main(sys.argv[1:])