### get parameters
args = commandArgs(trailingOnly=TRUE)
index_matrix_inputfile = args[1]
cREs_number_outfile = args[2]

### read index set matrix
read_matrix = function(inputfile){
	data_index_set = as.matrix(read.table(inputfile, header=T))
	rownames(data_index_set) = data_index_set[,1]
	data_index_set = data_index_set[,-1]
	class(data_index_set) = "numeric" 
	return(data_index_set)
}

###########
data_index_matrix = read_matrix(index_matrix_inputfile)

cRE_num = colSums(data_index_matrix) 

print('cREs number: ')
print(cRE_num)

pdf(cREs_number_outfile)
plot(cRE_num, xaxt='n', pch = 19, xlab = '')
axis(1, at=c(1:length(cRE_num)), labels=FALSE )

text(x=c(1:length(cRE_num)), par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]), labels=colnames(data_index_matrix) , srt=45, adj=1, xpd=TRUE)
dev.off()


# Rscript plot_index_set_region_hist.R celltype.index.sorted.txt cRE_celltype_number.pdf

