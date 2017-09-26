import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

################################################################################################
### read 2d array
def read2d_array(filename,dtype_used):
	import numpy as np
	data=open(filename,'r')
	data0=[]
	for records in data:
		tmp = [x.strip() for x in records.split('\t')]
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


index_set = read2d_array('celltype.index_set.sorted.txt', str)

index_set_data = np.log2(np.array(index_set[1:1000, 1:], dtype = float)+1)



my_cmap = mcolors.LinearSegmentedColormap.from_list(name='white_sth',colors=['white', 'black'],N=99)

plt.imshow(index_set_data, cmap=my_cmap, interpolation='nearest', vmin=0, vmax=15)
#plt.show()

plt.savefig('celltype.index_set.sorted.pdf', bbox_inches='tight', pad_inches=0, dpi=300)


