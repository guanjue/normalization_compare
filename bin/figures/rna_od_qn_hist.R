library(plotrix)
### get parameters
args = commandArgs(trailingOnly=TRUE)
qn_normed_inputfile = args[1]
od_sig_inputfile = args[2]
hist_outfile_qn = args[3]
hist_outfile_od_sig = args[4]

###########
### qn normalized matrix
qn_matrix = read.table(qn_normed_inputfile, header = F)
#print(summary(qn_matrix))

### qn normalized matrix
od_sig_matrix = read.table(od_sig_inputfile, header = F)
#print(summary(od_sig_matrix))

### plot histogram
pdf(hist_outfile_qn)
par(mfrow=c(2,2))
plot(density(log2(qn_matrix[,1]+0.01)), ylim=c(0,1), main='')
plot(density(log2(qn_matrix[,2]+0.01)), ylim=c(0,1), main='')
plot(density(log2(qn_matrix[,3]+0.01)), ylim=c(0,1), main='')
plot(density(log2(qn_matrix[,4]+0.01)), ylim=c(0,1), main='')
dev.off()

### plot histogram
pdf(hist_outfile_od_sig)
par(mfrow=c(2,2))
plot(density(log2(od_sig_matrix[,1]+0.01)), ylim=c(0,1), main='')
plot(density(log2(od_sig_matrix[,2]+0.01)), ylim=c(0,1), main='')
plot(density(log2(od_sig_matrix[,3]+0.01)), ylim=c(0,1), main='')
plot(density(log2(od_sig_matrix[,4]+0.01)), ylim=c(0,1), main='')
dev.off()

# Rscript rna_od_qn_hist.R all_rna_sampe_target_class_qn.txt all_rna_sampe_target_class.txt all_rna_sampe_target_class_qn_hist.png all_rna_sampe_target_class_hist.png

