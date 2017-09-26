import numpy as np
from scipy import stats

def tabular2composite(input_file_list,startcol):
	file_name_list=open(input_file_list,'r')

	sig_matrix_dict={}
	sig_matrix_rank_dict={}

	file_list=[]
	for file_0 in file_name_list:
		file_0=file_0.split()[0]
		print(file_0)
		file_list.append(file_0)
		data=open(file_0,'r')
		header=data.readline()
		sig_matrix=[]
		name_matrix=[]
		for sig in data:
			sig_tmp=[x.strip() for x in sig.split('\t')]
			name=sig_tmp[0:startcol-1]
			sig_matrix.append(sig_tmp[startcol-1:])
			name_matrix.append(name)

		sig_matrix=np.array(sig_matrix,dtype=float)
		sig_matrix_rank=stats.rankdata(sig_matrix, method='ordinal')
		sig_matrix_rank_dict[file_0]=sig_matrix_rank

		for r,s in zip(sig_matrix_rank,sig_matrix.flatten()):
			#print(r)
			if r in sig_matrix_dict:
				sig_matrix_dict[r].append(s)
			else:
				sig_matrix_dict[r]=[s]

	print('read data DONE!')
	quantile_norm={}
	j=0
	for i in sig_matrix_dict:
		j=j+1
		if j%100000==0:
			print(j)
		quantile_norm[i]=np.mean(sig_matrix_dict[i])

	print('Quantule Normalization DONE!')
	for file_1 in file_list:
		print(file_1)
		sig_matrix_qn_flatten=[]
		for i in sig_matrix_rank_dict[file_1]:
			sig_matrix_qn_flatten.append(quantile_norm[i])
		sig_matrix_qn_flatten=np.array(sig_matrix_qn_flatten)

		sig_matrix_qn=sig_matrix_qn_flatten.reshape(np.shape(sig_matrix))

		data_qn=open(file_1+'.qn.tabular','w')
		data_qn.write(header)
		for n,sig_qn in zip(name_matrix, sig_matrix_qn):
			for n in name:
				data_qn.write(n+'\t')
			for i in range(0,len(sig_qn)-1):
				data_qn.write(str(sig_qn[i])+'\t')
			data_qn.write(str(sig_qn[len(sig_qn)-1])+'\n')

		data_qn.close()
		data.close()

############################################################################
# python tabular_subtractnotag.py -i ref_pip_pkcenter_up1_1001_distsort_Htz1_i5006_BY4741_-_YPD_-_XO111_npf5017-160817_Pugh63322sacCer3_read1_combined.tabular -c 3 -o ref_pip_pkcenter_up1_1001_distsort_Htz1_read1_combined_bgsubtract.tabular -n ref_pip_pkcenter_up1_1001_distsort_masterNotag_20170213_read1_combined.tabular
import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:c:")
	except getopt.GetoptError:
		print 'python relocate_TSS.py -b bed_file -r readscount_file -o output_file -t searchTSS_win -w output_win'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'python relocate_TSS.py -b bed_file -r readscount_file -o output_file -t searchTSS_win -w output_win'
			sys.exit()
		elif opt=="-i":
			input_file_list=str(arg.strip())
		elif opt=="-c":
			startcol=int(arg.strip())

	tabular2composite(input_file_list,startcol)

if __name__=="__main__":
	main(sys.argv[1:])


