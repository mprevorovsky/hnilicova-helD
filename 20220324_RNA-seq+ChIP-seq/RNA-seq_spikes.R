chrom_sizes <- read.delim('chrom.sizes', header = FALSE, sep = '\t')
spike_reads <- read.delim('spike-in_read_counts_Msmegmatis', header = TRUE, sep = '\t', row.names = 1)
spike_reads <- spike_reads[c(-8:-11), ]
rownames(spike_reads) <- c('WT_EXP_1', 'helD_KO_EXP_3', 'WT_STAT_3', 'helD_KO_STAT_3', 'helD_KO_EXP_1', 'WT_STAT_1', 'helD_KO_STAT_1', 'WT_EXP_3')
barplot(spike_reads$other_mapped / spike_reads$spikes, names.arg = rownames(spike_reads), las = 3, cex.names = 0.6)

counts <- read.csv('counts_RNA-seq_Msmegmatis.csv', row.names = 1)
counts <- counts[, !colnames(counts) %in% c('WT_EXP_2', 'WT_STAT_2', 'helD_KO_STAT_2')]

plot(counts$WT_EXP_1, counts$WT_EXP_3, log = 'xy')
abline(a = 0, b = 1, col= 'blue')
plot(counts$WT_STAT_1, counts$WT_STAT_3, log = 'xy')
abline(a = 0, b = 1, col= 'blue')

plot(counts$helD_KO_EXP_1, counts$helD_KO_EXP_3, log = 'xy')
abline(a = 0, b = 1, col= 'blue')
plot(counts$helD_KO_STAT_1, counts$helD_KO_STAT_3, log = 'xy')
abline(a = 0, b = 1, col= 'blue')


counts.spike_norm <- counts[, -1]
for (i in colnames(counts.spike_norm)) {
  counts.spike_norm[, i] <- counts.spike_norm[, i] / spike_reads[i, 'spikes']
}

plot(counts.spike_norm$WT_EXP_1, counts.spike_norm$helD_KO_EXP_1, log = 'xy')
abline(a = 0, b = 1, col= 'blue')
plot(counts.spike_norm$WT_EXP_3, counts.spike_norm$helD_KO_EXP_3, log = 'xy')
abline(a = 0, b = 1, col= 'blue')

plot(counts.spike_norm$WT_STAT_1, counts.spike_norm$helD_KO_STAT_1, log = 'xy')
abline(a = 0, b = 1, col= 'blue')
plot(counts.spike_norm$WT_STAT_3, counts.spike_norm$helD_KO_STAT_3, log = 'xy')
abline(a = 0, b = 1, col= 'blue')

boxplot(log2(counts.spike_norm$helD_KO_EXP_1 / counts.spike_norm$WT_EXP_1),
        log2(counts.spike_norm$helD_KO_EXP_3 / counts.spike_norm$WT_EXP_3),
        log2(counts.spike_norm$helD_KO_STAT_1 / counts.spike_norm$WT_STAT_1),
        log2(counts.spike_norm$helD_KO_STAT_3 / counts.spike_norm$WT_STAT_3),
        ylim = c(-3, 3), names = c('EXP_1', 'EXP_3', 'STAT_1', 'STAT_3'), las = 2, 
        main = 'helD_KO / WT (spike-normalized)', ylab = 'log2(expression ratio)')
abline(h = 0, col = 'blue')

boxplot(log2(counts$helD_KO_EXP_1 / counts$WT_EXP_1),
        log2(counts$helD_KO_EXP_3 / counts$WT_EXP_3),
        log2(counts$helD_KO_STAT_1 / counts$WT_STAT_1),
        log2(counts$helD_KO_STAT_3 / counts$WT_STAT_3),
        ylim = c(-3, 3), names = c('EXP_1', 'EXP_3', 'STAT_1', 'STAT_3'), las = 2, 
        main = 'helD_KO / WT (raw)', ylab = 'log2(expression ratio)')
abline(h = 0, col = 'blue')


deseq.exp <- read.csv('DEG_RNA-seq_Msmegmatis/DESeq2results_helD-EXP_vs_WT-EXP.csv', row.names = 1)
plot(deseq.exp$baseMean, deseq.exp$log2FoldChange, log = 'x', pch = 20, 
     main = 'EXP', xlab = 'log(baseMean expression)', ylab = 'log2(foldChange expression)')
abline(h = c(1, -1), col = 'blue')
deseq.stat <- read.csv('DEG_RNA-seq_Msmegmatis/DESeq2results_helD-STAT_vs_WT-STAT.csv', row.names = 1)
plot(deseq.stat$baseMean, deseq.stat$log2FoldChange, log = 'x', pch = 20, 
     main = 'STAT', xlab = 'log(baseMean expression)', ylab = 'log2(foldChange expression)')
abline(h = c(1, -1), col = 'blue')

