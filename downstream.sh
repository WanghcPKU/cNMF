library(Seurat)
library(ggplot2)

iso.p8 <- readRDS("df.rds")
##iso.p8 已经运行过Seurat上游了
usage_file <- read.table("/lustre/user/taowlab/wanghc/work/lvwc/ytgs/zuhui/20240619/cNMF/isop8.fib_cNMF/isop8.fib_cNMF.usages.k_4.dt_0_1.consensus.txt",
sep='\t', row.names=1, header=TRUE)
head(usage_file)
spectra_score_file <- read.table("/lustre/user/taowlab/wanghc/work/lvwc/ytgs/zuhui/20240619/cNMF/isop8.fib_cNMF/isop8.fib_cNMF.gene_spectra_score.k_4.dt_0_1.txt",
sep='\t', row.names=1, header=TRUE)

usage_norm <- as.data.frame(t(apply(usage_file, 1, function(x) x / sum(x)))) ###标准化
head(usage_norm)

##合并矩阵
new_metadata <- merge(iso.p8@meta.data, usage_norm, by = "row.names", all.x = TRUE)
rownames(new_metadata) <- new_metadata$Row.names
iso.p8@meta.data <- new_metadata

##Vlnplot 确认所有的Program都有意义

usage_norm_t <- t(usage_norm)
data_long <- melt(usage_norm_t, id.vars = rownames(usage_norm_t))

pdf('./vlnplot.pdf') ###图片待美化
ggplot(data_long, aes(x = Var1, y = Values, fill = Var1)) + 
  geom_violin(trim = FALSE) + geom_jitter()+
  labs(title = "Violin Plot", x = "Cells", y = "Values") +
  theme_minimal()
dev.off()

###选择每个program里面的代表基因 这里选择top50的基因 选择后可以进行GO分析
get_top_colnames <- function(row) {
  print(row[1:5])
  top_indices <- order(row, decreasing = TRUE)[1:50]
  return(colnames(spectra_score_file)[top_indices])
}

top_colnames <- apply(spectra_score_file, 1, get_top_colnames)
top_colnames <- as.data.frame(top_colnames)

top_colnames

write.csv(top_colnames,"cNMF_marker.csv")

#### 可视化marker 与 program
top.matrix <- spectra_score_file[,unique(c(top_colnames[,1],top_colnames[,2],top_colnames[,3],top_colnames[,4]))]
markers <- unique(c(top_colnames[,1],top_colnames[,2],top_colnames[,3],top_colnames[,4]))
markers <- gsub("\\.", "-", markers)
iso.p8.scale.data <- iso.p8@assays$RNA@scale.data
head(iso.p8.scale.data)

top.matrix <- iso.p8.scale.data[markers,]

pdf('./test1.pdf')
pheatmap(top.matrix,scale="row",show_colnames = FALSE,show_rownames = FALSE,cluster_cols = FALSE,cluster_rows = FALSE)
dev.off()

### 可视化单细胞中不同program的表达情况
p <- FeaturePlot(iso.p8, features = colnames(usage_norm), combine=F)
pdf('./test2.pdf')
p
dev.off()
