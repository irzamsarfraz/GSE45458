

matrix <- readRDS("./data/matrix.rds")

sce <- SingleCellExperiment(list(counts = matrix))

gse45458 <- getGEO('GSE45458')

colData <- DataFrame(
  cell_type =  stringr::str_split(
    string = gse45458$GSE45458_series_matrix.txt.gz$characteristics_ch1, 
    pattern = "cell type: ", n = 2, simplify = TRUE)[, 2], 
  cytokine_treatment = stringr::str_split(
    string = gse45458$GSE45458_series_matrix.txt.gz$characteristics_ch1.1, 
    pattern = "cytokine treatment: ", n = 2, simplify = TRUE)[, 2],
  time_point = stringr::str_split(
    string = gse45458$GSE45458_series_matrix.txt.gz$characteristics_ch1.2, 
    pattern = "time point: ", n = 2, simplify = TRUE)[, 2],
  row.names = colnames(sce))

colData(sce) <- colData

# ILC2_Media_24_hour (4 replicates) vs. ILC2_IL-33_24_hour (4 replicates)
design <- cbind("ILC2_media_24_hour" = c(rep(1, 4), rep(0, 4)), "ILC2_IL_33_24_hour" = c(rep(0,4), rep(1,4)))
rownames(design) <- c(colnames(sce)[which(sce$cytokine_treatment == "media" & sce$time_point == "24 hours")], colnames(sce)[which(sce$cytokine_treatment == "IL-33" & sce$time_point == "24 hours")])
matrix <- assay(sce, "counts")[, rownames(design)]
fit <- lmFit(object = matrix, design = design)
cont.matrix <- makeContrasts(ILC2_Media_24_hour_vs_ILC2_IL_33_24_hour=ILC2_media_24_hour-ILC2_IL_33_24_hour, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
which(topTable(fit2, number = nrow(sce))$adj.P.Val < 0.05)
