## ----load, error=TRUE---------------------------------------------------------
trapnell.path <- "myoblast-differentiation_trapnell.rds"
if (!file.exists(trapnell.path)) {
  download.file("https://zenodo.org/record/1443566/files/real/gold/myoblast-differentiation_trapnell.rds?download=1",
                trapnell.path,
                mode="wb")
}
trapnell.dynverse <- readRDS(trapnell.path)

## ----preprocessing, error=TRUE------------------------------------------------
trapnell <- Seurat::CreateSeuratObject(counts=t(trapnell.dynverse$count),
                                       min.cells=3,
                                       min.features=200)

## ----visualization, error=TRUE, fig.asp=1, fig.cap="Figure 1. The result of PCA for the myoblast dataset"----
trapnell <- Seurat::FindVariableFeatures(trapnell, verbose=FALSE)
trapnell <- Seurat::ScaleData(trapnell, verbose=FALSE)
trapnell <- Seurat::RunPCA(trapnell, verbose=FALSE)
plot(Seurat::Embeddings(trapnell))

## ----estimate, error=TRUE-----------------------------------------------------
trapnell.fit <- treefit::treefit(trapnell)
trapnell.fit

## ----plot, error=TRUE, fig.cap="Figure 2. Summary of the Treefit analysis results"----
plot(trapnell.fit)

