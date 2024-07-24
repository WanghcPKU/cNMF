df <- readRDS("./df.rds")
df.count_matrix <- df@assays$RNA@counts
df.count_matrix_filtered <- df.count_matrix[rowSums(df.count_matrix) > 0, ] ## 过滤全为0

writeMM(df.count_matrix_filtered,"./cNMF/matrix.mtx") ###保存稀疏矩阵

barcodes <- colnames(df.count_matrix_filtered)
write.table(as.data.frame(barcodes), "./cNMF/barcodes.tsv",
           col.names = FALSE, row.names = FALSE, sep = "\t") 

gene_names <- rownames(df.count_matrix_filtered)
features <- data.frame("gene_id" = gene_names,"gene_name" = gene_names,type = "Gene Expression")
write.table(as.data.frame(features), sep = "\t", "./cNMF/genes.tsv",
           col.names = FALSE, row.names = FALSE)

