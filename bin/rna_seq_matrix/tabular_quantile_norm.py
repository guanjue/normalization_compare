import numpy as np
from scipy import stats

def tabular2composite(input_file,input_file_notag,startcol,output_file):
	data=open(input_file,'r')
	data_notag=open(input_file_notag,'r')
	data_qn=open(output_file,'w')
	data_bg_qn=open(output_file+'.bg.tabular','w')
	data_bg_sub_qn=open(output_file+'.bg_sub.tabular','w')

	header=data.readline()
	data_qn.write(header)
	data_bg_qn.write(data_notag.readline())
	data_bg_sub_qn.write(header)

	sig_matrix=[]
	bg_matrix=[]
	name_matrix=[]

	for sig,bg in zip(data,data_notag):
		sig_tmp=[x.strip() for x in sig.split('\t')]
		bg_tmp=[x.strip() for x in bg.split('\t')]

		name=sig_tmp[0:startcol-1]
		sig_matrix.append(sig_tmp[startcol-1:])
		bg_matrix.append(bg_tmp[startcol-1:])
		name_matrix.append(name)

	sig_matrix=np.array(sig_matrix,dtype=float)
	bg_matrix=np.array(bg_matrix,dtype=float)
	print(np.shape(sig_matrix))
	sig_matrix_rank=stats.rankdata(sig_matrix, method='ordinal')
	bg_matrix_rank=stats.rankdata(bg_matrix, method='ordinal')

	sig_matrix_dict={}
	bg_matrix_dict={}
	for r,s in zip(sig_matrix_rank,sig_matrix.flatten()):
		#print(r)
		sig_matrix_dict[r]=s
	for r,s in zip(bg_matrix_rank,bg_matrix.flatten()):
		bg_matrix_dict[r]=s

	#print((sig_matrix_dict)[3])
	quantile_norm={}
	for i in sig_matrix_rank:
		qn=(sig_matrix_dict[i]+bg_matrix_dict[i])/2
		quantile_norm[i]=qn

	sig_matrix_qn_flatten=[]
	for i in sig_matrix_rank:
		sig_matrix_qn_flatten.append(quantile_norm[i])
	sig_matrix_qn_flatten=np.array(sig_matrix_qn_flatten)

	sig_matrix_qn=sig_matrix_qn_flatten.reshape(np.shape(sig_matrix))


	bg_matrix_qn_flatten=[]
	for i in bg_matrix_rank:
		bg_matrix_qn_flatten.append(quantile_norm[i])
	bg_matrix_qn_flatten=np.array(bg_matrix_qn_flatten)

	bg_matrix_qn=bg_matrix_qn_flatten.reshape(np.shape(bg_matrix))


	for n,sig_qn in zip(name_matrix, sig_matrix_qn):
		for n in name:
			data_qn.write(n+'\t')
		for s in sig_qn:
			data_qn.write(str(s)+'\t')

		data_qn.write('\n')


	for n,sig_qn in zip(name_matrix, bg_matrix_qn):
		for n in name:
			data_bg_qn.write(n+'\t')
		for s in sig_qn:
			data_bg_qn.write(str(s)+'\t')

		data_bg_qn.write('\n')

	sig_matrix_qn_bg_sub=sig_matrix_qn-bg_matrix_qn
	for n,sig_qn in zip(name_matrix, sig_matrix_qn_bg_sub):
		for n in name:
			data_bg_sub_qn.write(n+'\t')
		for s in sig_qn:
			data_bg_sub_qn.write(str(s)+'\t')

		data_bg_sub_qn.write('\n')

	data_bg_sub_qn.close()
	data_qn.close()
	data.close()
	data_notag.close()
	data_bg_qn.close()

############################################################################
# python tabular_subtractnotag.py -i ref_pip_pkcenter_up1_1001_distsort_Htz1_i5006_BY4741_-_YPD_-_XO111_npf5017-160817_Pugh63322sacCer3_read1_combined.tabular -c 3 -o ref_pip_pkcenter_up1_1001_distsort_Htz1_read1_combined_bgsubtract.tabular -n ref_pip_pkcenter_up1_1001_distsort_masterNotag_20170213_read1_combined.tabular
import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:n:c:o:")
	except getopt.GetoptError:
		print 'python relocate_TSS.py -b bed_file -r readscount_file -o output_file -t searchTSS_win -w output_win'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'python relocate_TSS.py -b bed_file -r readscount_file -o output_file -t searchTSS_win -w output_win'
			sys.exit()
		elif opt=="-i":
			input_file=str(arg.strip())
		elif opt=="-n":
			input_file_notag=str(arg.strip())
		elif opt=="-c":
			startcol=int(arg.strip())
		elif opt=="-o":
			output_file=str(arg.strip())

	tabular2composite(input_file,input_file_notag,startcol,output_file)

if __name__=="__main__":
	main(sys.argv[1:])


