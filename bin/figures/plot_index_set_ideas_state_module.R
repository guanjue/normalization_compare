library(pheatmap)
library(Matrix)

### get parameters
args = commandArgs(trailingOnly=TRUE)
index_set_inputfile = args[1]
index_set_se_inputfile = args[2]

index_inputfile = args[3]
index_set_inputfile_binary = args[4]

index_set_all_heatmap = args[5]
index_set_thresh_heatmap = args[6]

index_set_se_all_heatmap = args[7]
index_set_se_thresh_heatmap = args[8]

index_all_heatmap = args[9]

color = args[10]
threshold = as.numeric(args[11])
upper_lim_quantile = as.numeric(args[12])
tranform = args[13]


### read index set signal matrix
read_index_set_sig_matrix = function(inputfile){
	data_index_set = as.matrix(read.table(inputfile, header=T))
	rownames(data_index_set) = data_index_set[,1]
	data_index_set = data_index_set[,-1]
	class(data_index_set) = "numeric" 
	return(data_index_set)
}

### read index set matrix
read_matrix = function(inputfile){
	data_index_set = as.matrix(read.table(inputfile, header=T))
	rownames(data_index_set) = data_index_set[,1]
	data_index_set = data_index_set[,-1]
	class(data_index_set) = "numeric" 
	return(data_index_set)
}

###########
### read index set signal matrix
data_index_set_sig = read_index_set_sig_matrix(index_set_inputfile)
### read index set shannon entropy matrix
data_index_set_se = read_index_set_sig_matrix(index_set_se_inputfile)

### read index matrix
data_index = read_matrix(index_inputfile)
### read index set binary pattern matrix
data_index_set_01 = read_matrix(index_set_inputfile_binary)

### set heatmap colors
rgb_col_num=read.table('input_data/state_color.txt',header=F)
rgb_col_num=rbind(rgb_col_num, c(255,255,255))
print(rgb_col_num)
rgb_col=apply(rgb_col_num,1,function(x) rgb(x[1],x[2],x[3],max=255))

my_colorbar=colorRampPalette(rgb_col)(n = 35)
col_breaks = c(
        seq(0.1, 1,length=2),
        seq(1.1, 2,length=2),
        seq(2.1, 3,length=2),
        seq(3.1, 4,length=2),
        seq(4.1, 5,length=2),
        seq(5.1, 6,length=2),
        seq(6.1, 7,length=2),
        seq(7.1, 8,length=2),
        seq(8.1, 9,length=2),
        seq(9.1, 10,length=2),
        seq(10.1, 11,length=2),
        seq(11.1, 12,length=2),
        seq(12.1, 13,length=2),
        seq(13.1, 14,length=2),
        seq(14.1, 15,length=2),
        seq(15.1, 16,length=2),
        seq(16.1, 17,length=2),
        seq(17.1, 18,length=2)
)

### plot index set heatmap
pheatmap(data_index_set_sig+1, color=my_colorbar, breaks=col_breaks, cluster_cols = FALSE,cluster_rows=FALSE,annotation_names_row=FALSE,annotation_names_col=TRUE,show_rownames=FALSE,show_colnames=TRUE, filename = index_set_all_heatmap)
### plot index set heatmap (filter by DNA region number threshold)
data_index_set_filtered = data_index_set_sig[apply(data_index_set_01, 1, max) >= threshold, ]
pheatmap(data_index_set_filtered+1, color=my_colorbar, breaks=col_breaks, cluster_cols = FALSE,cluster_rows=FALSE,annotation_names_row=FALSE,annotation_names_col=TRUE,show_rownames=FALSE,show_colnames=TRUE, filename = index_set_thresh_heatmap)

### plot index set Shannon Entropy heatmap
pheatmap(data_index_set_se, color=colorRampPalette(c(rgb(1,1,1,1),rgb(1,1,1,0)), alpha = TRUE)(n = 128), cluster_cols = FALSE,cluster_rows=FALSE,annotation_names_row=FALSE,annotation_names_col=TRUE,show_rownames=FALSE,show_colnames=TRUE, filename = index_set_se_all_heatmap)
### plot index set Shannon Entropy heatmap (filter by DNA region number threshold)
data_index_set_se_filtered = data_index_set_se[apply(data_index_set_01, 1, max) >= threshold, ]
pheatmap(data_index_set_se_filtered, color=colorRampPalette(c(rgb(1,1,1,1),rgb(1,1,1,0)), alpha = TRUE)(n = 128), cluster_cols = FALSE,cluster_rows=FALSE,annotation_names_row=FALSE,annotation_names_col=TRUE,show_rownames=FALSE,show_colnames=TRUE, filename = index_set_se_thresh_heatmap)

###########
### plot index heatmap
pheatmap(data_index+1, color=my_colorbar, breaks=col_breaks, cluster_cols = FALSE, cluster_rows=FALSE,annotation_names_row=FALSE,annotation_names_col=TRUE,show_rownames=FALSE,show_colnames=TRUE, filename = index_all_heatmap)
###########

# Rscript plot_index_set_module.R celltype.index_set.sorted.txt celltype.binary_pattern.sorted.txt black 200 index_set_all.pdf index_set_thresh.pdf index.png