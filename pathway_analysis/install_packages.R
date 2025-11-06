# Install CRAN packages
install.packages(c("shiny", "DT", "ggplot2", "dplyr"))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("clusterProfiler",
                       "org.Hs.eg.db",
                       "org.Mm.eg.db",
                       "enrichplot"))
