#Run Rv4.1.0

######### --Libraries--#########
options(mc.cores = 8)
library(limma)
library(edgeR)
library(stats)
library(factoextra)
library(umap)
library(ggplot2)
library(ExpressionNormalizationWorkflow)
library(stringr)
library(pvca)
library(mclust)
library(fossil)
library(amap)
library(dendextend)
library(data.table)
library(rrcov)
library(svglite)
library("Hmisc")
library(corrplot)
source("./heatmap3LW_function.r")

# TITLE: Sacha 01 CyCMV flu vaccine
# AUTHOR: LEANNE WHITMORE

######### --FILES/OBJECTS--#########
rhesus2human <- read.csv(
    file = "rhesus2human.csv",
    header = TRUE,
    stringsAsFactors = FALSE
)

######### --FUNCTIONS--#########

theme_minimal_LW <- function(base_size = 14, base_family = "arial") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size = base_size, base_family = base_family)
    + theme(
            plot.title = element_text(
                face = "bold",
                size = rel(1.2), hjust = 0.5
            ),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold", size = rel(1)),
            axis.title.y = element_text(angle = 90, vjust = 2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size = unit(0.4, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(size = 10, face = "bold"),
            plot.margin = unit(c(10, 5, 5, 5), "mm"),
            strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
            strip.text = element_text(face = "bold")
        ))
}

theme_Publication <- function(base_size = 14, base_family = "arial") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size = base_size, base_family = base_family)
    + theme(
            plot.title = element_text(
                face = "bold",
                size = rel(1.2), hjust = 0.5
            ),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold", size = rel(1)),
            axis.title.y = element_text(angle = 90, vjust = 2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour = "#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size = unit(0.6, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face = "italic"),
            plot.margin = unit(c(10, 5, 5, 5), "mm"),
            strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
            strip.text = element_text(face = "bold")
        ))
}

scale_fill_Publication <- function(..) {
    library(scales)
    discrete_scale("fill", "Publication", manual_pal(values = c("#386cb0", "#fdb462", "#7fc97f", "#ef3b2c", "#662506", "#a6cee3", "#fb9a99", "#984ea3", "#ffff33")), ..)
}

scale_colour_Publication <- function(..) {
    library(scales)
    discrete_scale("colour", "Publication", manual_pal(values = c("#386cb0", "#fdb462", "#7fc97f", "#ef3b2c", "#662506", "#a6cee3", "#fb9a99", "#984ea3", "#ffff33")), ..)
}

generate_folder <- function(foldername) {
    workDir <- getwd()
    subDir <- foldername
    results_path <- file.path(workDir, subDir)
    if (file.exists(subDir)) {
    } else {
        dir.create(results_path)
    }
    return(results_path)
}

vizualize_DE_genes_bp <- function(results, plot_file) {
    print("STATUS: Generating bar plot of number of DE genes...")
    results_t <- t(summary(results))
    results_t <- results_t[, -2]

    for (i in 1:(length(row.names(results_t)))) {
        results_t[i, 1] <- results_t[i, 1] * -1
    }

    DE <- as.data.frame(results_t)
    DE <- setnames(DE,
        old = c("Var1", "Var2", "Freq"),
        new = c("Time_Point", "group", "DE_genes")
    )

    # Create plot
    ggplot(DE, aes(
        x = Time_Point, y = DE_genes, fill = group,
        label = DE$DE_genes
    )) +
        geom_bar(stat = "identity", position = "identity") +
        # geom_text(size = 5, position = position_stack(vjust = 0) )+
        # theme_light() +
        theme_minimal() +
        scale_fill_manual(values = c("#0808c4", "#da9618")) +
        # xlab("Time point")
        ylab("Number of Differentially Expressed Genes") +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
            axis.text.y = element_text(size = 15)
        )
    ggsave(plot_file, width = 6, height = 4, units="in", dpi = 300)
}

generate_boxplots_voom <- function(data, labels, filename, figres, maintitle, ylabtitle) {
    png(filename, width = 10, height = 8, units = "in", res = figres)
    # par(mar=c(1,1,1,1))
    minvalue <- min(data)
    maxvalue <- max(data)
    boxplot(data,
        labels = labels, ylim = c(minvalue - 1, maxvalue + 1),
        ylab = ylabtitle, main = maintitle, cex.axis = .6, las = 2,
        frame = FALSE
    )
    dev.off()
}

generate_density_plot <- function(data, labels, filename, figres) {
    png(filename, res = figres)
    par(xpd = TRUE)
    if (length(labels) > 10) {
        plotDensities(data, legend = FALSE)
    } else {
        plotDensities(data,
            legend = "topright",
            inset = c(-0.2, 0), levels(labels)
        )
    }
    dev.off()
}

normalize_data <- function(CM2, targetfile, group, control = FALSE) {

    # order target and count matrix so they are the same (THIS IS IMPORTANT)
    CM2 <- CM2[, rownames(targetfile)]

    # CHECK IF ORDER IS THE SAME
    if (all.equal(colnames(CM2), rownames(targetfile)) != TRUE) {
        print("MASSIVE WARNING: RESULTS WILL BE WRONG IF THIS IS NOT EQUAL!!!!!!!!")
        print(rownames(targetfile))
        print(colnames(CM2))
    }

    # normalize
    CM2 <- DGEList(counts = CM2)
    CM2 <- calcNormFactors(CM2, method = "TMM") # TMM normalization
    png(file.path(norm_results, "mean_variance_norm.png"))
    Pi.CPM <- voom(counts = CM2, normalize.method = "none", plot = T, span = 0.1)
    dev.off()
    write.csv(Pi.CPM$E, file.path(norm_results, paste0("1.norm_matrix_", group, ".csv")))

    sig_HGNC <- merge(rhesus2human, Pi.CPM$E,
        by.x = "Gene.stable.ID",
        by.y = "row.names",
        all.X = T, all.Y = T
    )

    sig_HGNC <- sig_HGNC[, !(names(sig_HGNC) %in% c("Gene.stable.ID"))]
    sig_HGNC <- avereps(sig_HGNC,
        ID = sig_HGNC$HGNC.symbol
    )
    rownames(sig_HGNC) <- sig_HGNC[, "HGNC.symbol"]
    sig_HGNC <- sig_HGNC[, !(colnames(sig_HGNC) %in% c("HGNC.symbol"))]
    sig_HGNC <- as.matrix(data.frame(sig_HGNC))
    write.csv(sig_HGNC, file.path(norm_results, paste0("1.norm_matrix_HGNC_", group, ".csv")), quote = FALSE)
    return(Pi.CPM)
}

pca_fun <- function(exprs, labels, results_path,
                    base_file_name, target_columns,
                    figres = 100, size = 1, pca=FALSE, legend="right") {
    # Run PCA/SVD reduction
    if (isFALSE(pca)) {
        pca <- prcomp(t(exprs))
    }
    E <- get_eig(pca)
    cx <- sweep(t(exprs), 2, colMeans(t(exprs)), "-")
    sv <- svd(cx)


    vizualize_pca(
        file.path(results_path, paste0("svd_", base_file_name)),
        sv$u, labels[, target_columns[1]],
        labels[, target_columns[2]], figres, E, size, legend
    )
    vizualize_pca(
        file.path(results_path, paste0("pca_", base_file_name)),
        pca$x, labels[, target_columns[1]],
        labels[, target_columns[2]],
        figres, E, size, legend
    )
    vizualize_scree_plot(
        file.path(
            results_path,
            paste0("scree_", base_file_name)
        ), pca, figres
    )

    loadingscores <- as.data.frame(pca$rotation)
    is_pc1_0 <- loadingscores$PC1 > 0
    is_pc2_0 <- loadingscores$PC2 > 0

    loadingscores <- loadingscores[is_pc1_0, ]
    loadingscores <- loadingscores[with(loadingscores, order(-PC1)), ]
    save_loading_scores(
        file.path(results_path, paste0("loadingscores_pc1", base_file_name, ".txt")),
        loadingscores["PC1"], figres
    )

    loadingscores <- as.data.frame(pca$rotation)
    loadingscores <- loadingscores[is_pc2_0, ]
    loadingscores <- loadingscores[with(loadingscores, order(-PC2)), ]
    save_loading_scores(
        file.path(results_path, paste0("loadingscores_pc2", base_file_name, ".txt")),
        loadingscores["PC2"], figres
    )
    return(pca)
}

umap_fun <- function(exprs, labels, results_path,
                     base_file_name, target_columns,
                     figres = 100, size = 1, UMAP=FALSE, legend="right") {
    # Runs default paramaters of umap
    if (isFALSE(UMAP)) {
        UMAP <- umap(t(exprs))
    }
    vizualize_umap(
        file.path(results_path, paste0("umap_", base_file_name)),
        UMAP$layout, labels[, target_columns[1]],
        labels[, target_columns[2]], figres, size, legend
    )

    return(UMAP)
}

vizualize_umap <- function(plot_file, U, class1, class2, figres, size, legend) {
    # Vizualize umap reduction
    library(Polychrome)
    minx <- min(U[, 1])
    maxx <- max(U[, 1])
    miny <- min(U[, 2])
    maxy <- max(U[, 2])
    if (length(levels(factor(class2))) <= 3) {
        if (length(levels(factor(class1))) <= 6) {
            qplot(U[, 1], U[, 2], shape = factor(paste(class1)), color = factor(class2), size = I(size)) +
                theme_Publication() + theme(legend.title = element_blank()) +
                xlab("UMAP 1") +
                ylab("UMAP 2") +
                xlim(minx, maxx) + ylim(miny, maxy) +
                scale_color_manual(values = c("Protected" = "red", "NonProtected" = "black", "Deceased" = "gray")) +
                scale_fill_manual(values =c("Protected" = "red", "NonProtected" = "black", "Deceased" = "gray")) +
                theme(legend.position = legend)
        } else {
            qplot(U[, 1], U[, 2], shape = factor(paste(class1)), color = factor(class2), size = I(size)) +
                theme_Publication() + theme(legend.title = element_blank()) +
                xlab("UMAP 1") +
                ylab("UMAP 2") +
                xlim(minx, maxx) + ylim(miny, maxy) +
                scale_color_manual(values = c("Protected" = "red", "NonProtected" = "black", "Deceased" = "gray")) +
                scale_fill_manual(values =c("Protected" = "red", "NonProtected" = "black", "Deceased" = "gray")) +
                theme(legend.position = legend) +
                scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
        }
    } else {
        P36 <- createPalette(length(levels(factor(class2))), c("#ff0000", "#00ff00", "#0000ff"))
        if (length(levels(factor(class1))) <= 6) {
            qplot(U[, 1], U[, 2], shape = factor(paste(class1)), color = factor(class2), size = I(size)) +
                theme_Publication() + theme(legend.title = element_blank()) +
                scale_color_manual(values = as.character(P36)) +
                scale_fill_manual(values = as.character(P36)) +
                xlab("UMAP 1") +
                ylab("UMAP 2") +
                xlim(minx, maxx) + ylim(miny, maxy) +
                theme(legend.position = legend)
        } else if (length(levels(factor(class1))) > 6) {
            qplot(U[, 1], U[, 2], shape = factor(paste(class1)), color = factor(class2), size = I(size)) +
                theme_Publication() + theme(legend.title = element_blank()) +
                xlab("UMAP 1") +
                ylab("UMAP 2") +
                scale_color_manual(values = as.character(P36)) +
                scale_fill_manual(values = as.character(P36)) +
                xlim(minx, maxx) + ylim(miny, maxy) +
                theme(legend.position = legend) +
                scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
        }
    }
    ggsave(plot_file, width = 6, height = 4, units = "in", dpi = figres)
}

vizualize_pca <- function(plot_file, PCA, class1, class2, figres, E, size, legend) {
    # Vizualize PCA  results
    library(Polychrome)
    minx <- min(PCA[, 1])
    maxx <- max(PCA[, 1])
    miny <- min(PCA[, 2])
    maxy <- max(PCA[, 2])
    if (length(levels(factor(class2))) <= 3) {
        if (length(levels(factor(class1))) <= 6) {
            qplot(PCA[, 1], PCA[, 2], color = factor(class2), shape = factor(class1), size = I(size)) +
                theme_Publication() +
                theme(legend.title = element_blank()) +
                scale_color_manual(values = c("Protected" = "red", "NonProtected" = "black", "Deceased" = "gray")) +
                scale_fill_manual(values =c("Protected" = "red", "NonProtected" = "black", "Deceased" = "gray")) +
                xlab(paste0("PC1 ", round(E$variance.percent[1], digits = 2), "%")) +
                ylab(paste0("PC2 ", round(E$variance.percent[2], digits = 2), "%")) +
                theme(legend.position = legend)
        } else {
            qplot(PCA[, 1], PCA[, 2], color = factor(class2), shape = factor(class1), size = I(size)) +
                theme_Publication() +
                theme(legend.title = element_blank()) +
                scale_color_manual(values = c("IFN" = "pink", "CTL" = "black")) +
                scale_fill_manual(values = c("IFN" = "pink", "CTL" = "black")) +
                xlab(paste0("PC1 ", round(E$variance.percent[1], digits = 2), "%")) +
                ylab(paste0("PC2 ", round(E$variance.percent[2], digits = 2), "%")) +
                theme(legend.position = legend) +
                scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
        }
    } else {
        P36 <- createPalette(length(levels(factor(class2))), c("#ff0000", "#00ff00", "#0000ff"))
        if (length(levels(factor(class1))) <= 6) {
            qplot(PCA[, 1], PCA[, 2], color = factor(class2), shape = factor(class1), size = I(size)) +
                theme_Publication() +
                theme(legend.title = element_blank()) +
                xlab(paste0("PC1 ", round(E$variance.percent[1], digits = 2), "%")) +
                ylab(paste0("PC2 ", round(E$variance.percent[2], digits = 2), "%")) +
                theme(legend.position = legend) +
                scale_color_manual(values = as.character(P36)) +
                scale_fill_manual(values = as.character(P36))
        } else if (length(levels(factor(class1))) > 6) {
            qplot(PCA[, 1], PCA[, 2], color = factor(class2), shape = factor(class1), size = I(size)) +
                theme_Publication() +
                theme(legend.title = element_blank()) +
                xlab(paste0("PC1 ", round(E$variance.percent[1], digits = 2), "%")) +
                ylab(paste0("PC2 ", round(E$variance.percent[2], digits = 2), "%")) +
                theme(legend.position = legend) +
                scale_color_manual(values = as.character(P36)) +
                scale_fill_manual(values = as.character(P36)) +
                scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
        }
    }
    ggsave(plot_file, width = 6, height = 4, units = "in", dpi = 300)
}

vizualize_scree_plot <- function(plot_file, PCA, figres) {
    # Vizualize principle component variation results
    scree.plot <- fviz_eig(PCA, addlabels = TRUE, hjust = -0.3)
    png(plot_file, width = 7, height = 6, units = "in", res = figres)
    print(scree.plot)
    dev.off()
}

save_loading_scores <- function(write_file, df, figres) {
    # Save list of genes that have a positive effect on variation of principle
    # component 1 and 2 sorted from most influential
    write.table(df, file = write_file)
}

filter_read_counts_mean <- function(cm, filter_cutoff) {
    # Filter value was calculated by:
    #### Filters by row means usually set at 10 reads per gene across all samples

    A <- rowMeans(cm)
    isexpr <- A >= filter_cutoff
    cmfl <- cm[isexpr, ]
    return(cmfl)
}

rename_samples <- function(samples) {
    newsampleIDs <- c()
    for (i in samples) {
        i <- str_remove(i, "_RNA\\d+_Lib\\d+\\S*$")
        i <- str_replace_all(i, "-", "_")
        newsampleIDs <- c(newsampleIDs, i)
    }
    return(newsampleIDs)
}

generate_design_matrix <- function(normmatrix, target, groups) {
    vax <- factor(target[, "Vaccine_Treatment"])
    ti <- factor(target[, "TimePoint"])
    ID <- factor(target[, "Animal_ID"])
    tr <- factor(target[, "Animal_Outcome"])
    sex <- factor(target[, "Sex"])
    mm <- model.matrix(~ 0 + vax:ti:tr)

    rownames(mm) <- colnames(normmatrix)
    colnames(mm) <- make.names(colnames(mm))
    mm <- mm[, colnames(mm)[order(tolower(colnames(mm[, ])))]]
    mm <- mm[, colSums(mm) > 0]

    excludeAll <- nonEstimable(mm)
    if (length(excludeAll) > 0) {
        message("WARNING: These samples are nonEstimatable, design matrix ", excludeAll)
    }

    if ("ti" %in% excludeAll) {
        return("interactions term non estimable")
    }
    mm <- mm[, !colnames(mm) %in% excludeAll]
    if (!is.fullrank(mm)) {
        return("not full rank")
    }

    return(mm)
}

generate_clusterdendograms <- function(hc, plotfilename1, adjvalue, labelsize = 0.7) {
    counter <- 0
    labelsf <- c()
    colors <- c()
    dend <- as.dendrogram(hc)
    dend_labels <- labels(dend)
    P36 <- createPalette(length(levels(factor(dend_labels))), c("#ff0000", "#00ff00", "#0000ff"))
    names(P36) <- unique(dend_labels)
    for (i in dend_labels) {
        labelsf <- c(labelsf, i)
        colors <- c(colors, P36[[i]])
    }
    labels_colors(dend) <- colors
    labels_cex(dend) <- labelsize
    png(plotfilename1,
        units = "in", # bg = "transparent",
        width = 14.5, height = 5, res = 300
    )
    par(mar = c(6, 3, 2, 0.5), xpd = TRUE)
    plot(dend, xlab = "", main = "")
    mtext(paste0("Adj Rand index ", round(adjvalue, 3)))
    dev.off()
}


######### --Load in data--#########

# --Read in target files
message("STATUS: Load tables")
cm <- read.table("../count_matrix.txt", header = TRUE, sep = "\t", row.names = 1, as.is = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
target <- read.csv("./targetfile.csv", sep = ",", row.names = 1, as.is = TRUE, check.names = FALSE, header = TRUE, stringsAsFactors = FALSE)

target$Animal_Outcome <- str_replace_all(target$Animal_Outcome, "Non protected", "NonProtected")
# --Rename sample names in count matrix
newsampleIDs <- rename_samples(colnames(cm))
colnames(cm) <- newsampleIDs

# # remove animal that is deceased from analysis
# target <- target[target$Animal_Outcome != "Deceased", ]
# cm <- cm[, rownames(target)]
if (all.equal(colnames(cm), rownames(target)) != TRUE) {
    print("MASSIVE WARNING: RESULTS WILL BE WRONG IF THIS IS NOT EQUAL!!!!!!!!")
    cm <- cm[, rownames(target)]
    if (all.equal(colnames(cm), rownames(target)) == TRUE) {
        print("ISSUE HAS BEEN FIXED")
    } else {
        print("ISSUE IS NOT FIXED, PLEASE LOOK AT MANUALLY")
    }
}

# --Write new renamed files to file
count_results <- "1.count_data/"
generate_folder(count_results)
write.csv(cm, file.path(count_results, "count_matrix_renamed.csv"))
write.csv(target, file.path(count_results, "target_renamed.csv"))

# --generate figure of all counts
generate_density_plot(
    cm, rownames(target), file.path(count_results, "de_intensities_raw_counts.png"), 100
)

# -- Normalize data
norm_results <- "1.norm_data/"
generate_folder(norm_results)
Pi.CPM <- normalize_data(cm, target, "all")
generate_boxplots_voom(Pi.CPM$E, target$GaleID,
    file.path(norm_results, "boxplot_vnorm_all.png"),
    100,
    maintitle = "Normalized count matrix",
    ylabtitle = "voom normalized expression"
)

#-- filter out genes from each group that are below mean count of 5 across samples 
#Iteratively adjusted thresholds and decided that 5 was the best cutoff to get rid of the 
# mean-variance hook shown in the initial mean-variance plot in 1.norm_results/ (iterative results not saved just ran in R)
cmfl_counts <- filter_read_counts_mean(cm, 3)
write.csv(cmfl_counts, file.path(count_results, "count_matrix_renamed_fl.csv"))

norm_results <- "1.norm_data_fl/"
generate_folder(norm_results)
Pi.CPM <- normalize_data(cmfl_counts, target, "all")
generate_boxplots_voom(Pi.CPM$E, target$GaleID,
    file.path(norm_results, "boxplot_vnorm_all.png"),
    100,
    maintitle = "Normalized count matrix",
    ylabtitle = "voom normalized expression"
)

saveRDS(Pi.CPM, file.path(norm_results,"normobject.rds"))
# Pi.CPM <- readRDS(file.path(norm_results, "normobject.rds"))
cibersort <- Pi.CPM$E
cibersort <- merge(rhesus2human, cibersort,
    by.x = "Gene.stable.ID",
    by.y = "row.names",
    all.X = T, all.Y = T
)
cibersort <- cibersort[, !(names(cibersort) %in% c("Gene.stable.ID"))]
cibersort <- avereps(cibersort,
  ID = cibersort$HGNC.symbol
)
rownames(cibersort) <- cibersort[, "HGNC.symbol"]
cibersort <- cibersort[, !(colnames(cibersort) %in% c("HGNC.symbol"))]
cibersort <- as.matrix(data.frame(cibersort))
# colnames(cibersort) <- paste(target$Patient, target$Visit.Number, sep="_") 
class(cibersort) <- "numeric"

write.table(cibersort, file.path(norm_results, "1.norm_HGNC_cibersort.txt"), sep="\t")

######### --Feature Reduction--#########
message("STATUS: Run Feature reduction")
feature_results <- "1.feature_red"
generate_folder(feature_results)

pca <- pca_fun(
    Pi.CPM$E, target,
    feature_results, "_norm.png",
    c("Sex", "Animal_Outcome"), 300, 3
)

umap <- umap_fun(
    Pi.CPM$E, target,
    feature_results, "_norm.png",
     c("Sex", "Animal_Outcome"), 300, 3
)

umap <- umap_fun(
    Pi.CPM$E, target,
    feature_results, "_AnimalIDnorm.png",
    c("Animal_ID", "Animal_Outcome"), 300, 3, UMAP=umap
)
umap <- umap_fun(
    Pi.CPM$E, target,
    feature_results, "_AnimalIDTimePointnorm.png",
    c("Animal_ID", "TimePoint"), 300, 3,
    UMAP = umap
)

######### --Hierarchal clustering--#########

message("STATUS: Running hierarchal clustering...")
hi_cluster_results <- "1.hi_cluster_results"
generate_folder(hi_cluster_results)
d <- Pi.CPM$E
colnames(d) <- target$Animal_ID
hc <- hcluster(t(d), method = "pearson", link = "average")
x <- cutree(hc, k = length(unique(target$Animal_ID)))
m <- adj.rand.index(as.integer(x), as.integer(factor(as.character(target$Animal_ID))))
generate_clusterdendograms(hc,
    paste0(hi_cluster_results, "/dendogram_Patient.png"), m, labelsize = 1)

colnames(d) <- paste(target$Animal_ID, target$TimePoint, sep = "_")
hc <- hcluster(t(d), method = "pearson", link = "average")
x <- cutree(hc, k = length(unique(target$Animal_ID)))
m <- adj.rand.index(as.integer(x), as.integer(factor(as.character(target$Animal_ID))))
generate_clusterdendograms(hc,
    paste0(hi_cluster_results, "/dendogram_Patienttime.png"), m,
    labelsize = 1
)

######### --de analysis--#########
deresults_path <- "1.de"
generate_folder(deresults_path)

mm_all <- generate_design_matrix(Pi.CPM$E, target, groups = c("all"))
dupcor <- duplicateCorrelation(Pi.CPM, mm_all, block = target$Animal_ID)
dupcor$consensus.correlation
Pi.lmfit <- lmFit(Pi.CPM, design = mm_all, block = target$Animal_ID, correlation = dupcor$consensus)

contrastsmatrix <- c(
    "vax68.1.CyCMV.flu.tiD3.trProtected -vax68.1.CyCMV.flu.tiD0.trProtected",
    "vax68.1.CyCMV.flu.tiD14.trProtected -vax68.1.CyCMV.flu.tiD0.trProtected",
    "vax68.1.CyCMV.flu.tiD105.trProtected -vax68.1.CyCMV.flu.tiD0.trProtected",
    "vax68.1.CyCMV.flu.tiD108.trProtected -vax68.1.CyCMV.flu.tiD0.trProtected",
    "vax68.1.CyCMV.flu.tiD119.trProtected -vax68.1.CyCMV.flu.tiD0.trProtected",
    "vax68.1.CyCMV.flu.tiD3.trNonProtected -vax68.1.CyCMV.flu.tiD0.trNonProtected",
    "vax68.1.CyCMV.flu.tiD14.trNonProtected -vax68.1.CyCMV.flu.tiD0.trNonProtected",
    "vax68.1.CyCMV.flu.tiD105.trNonProtected -vax68.1.CyCMV.flu.tiD0.trNonProtected",
    "vax68.1.CyCMV.flu.tiD108.trNonProtected -vax68.1.CyCMV.flu.tiD0.trNonProtected",
    "vax68.1.CyCMV.flu.tiD119.trNonProtected -vax68.1.CyCMV.flu.tiD0.trNonProtected",

    "vaxfull.length.CyCMV.flu.tiD3.trProtected -vaxfull.length.CyCMV.flu.tiD0.trProtected",
    "vaxfull.length.CyCMV.flu.tiD14.trProtected -vaxfull.length.CyCMV.flu.tiD0.trProtected",
    "vaxfull.length.CyCMV.flu.tiD105.trProtected -vaxfull.length.CyCMV.flu.tiD0.trProtected",
    "vaxfull.length.CyCMV.flu.tiD108.trProtected -vaxfull.length.CyCMV.flu.tiD0.trProtected",
    "vaxfull.length.CyCMV.flu.tiD119.trProtected -vaxfull.length.CyCMV.flu.tiD0.trProtected",
    "vaxfull.length.CyCMV.flu.tiD3.trNonProtected -vaxfull.length.CyCMV.flu.tiD0.trNonProtected",
    "vaxfull.length.CyCMV.flu.tiD14.trNonProtected -vaxfull.length.CyCMV.flu.tiD0.trNonProtected",
    "vaxfull.length.CyCMV.flu.tiD105.trNonProtected -vaxfull.length.CyCMV.flu.tiD0.trNonProtected",
    "vaxfull.length.CyCMV.flu.tiD108.trNonProtected -vaxfull.length.CyCMV.flu.tiD0.trNonProtected",
    "vaxfull.length.CyCMV.flu.tiD119.trNonProtected -vaxfull.length.CyCMV.flu.tiD0.trNonProtected"
)
contr <- makeContrasts(contrasts = contrastsmatrix, levels = mm_all)

fit <- contrasts.fit(Pi.lmfit, contr) # This results in a 3-dimensional array, which is not entirely analogous to Pi.contrasts from above. The next two lines extract the coefficients and p-value matrices from the array, which is all you need to use the decideTests() function.
fit2 <- eBayes(fit)
results <- decideTests(fit2, lfc = (.58), method = "separate", adjust.method = "BH", p.value = 0.05)
summary(results)

dataMatrixde <- fit2$coefficients
sigMask <- dataMatrixde * (results**2) # 1 if significant, 0 otherwise
ExpressMatrixde <- subset(dataMatrixde, rowSums(sigMask) != 0)

# filter for significant genes - up/down regulated
sigMask <- subset(sigMask, rowSums(sigMask) != 0)
# Pi.contrasts$genes <- data.frame(ID_REF=rownames(Pi.contrasts))

write.csv(ExpressMatrixde, file = file.path(deresults_path, "expression_matrix_de.csv"), quote = F)
write.csv(results, file = file.path(deresults_path, "results_de.csv"), quote = F)
write.csv(dataMatrixde, file = file.path(deresults_path, "full_expression_matrix_de.csv"), quote = F)

ExpressMatrixde_HGNC <- merge(rhesus2human, ExpressMatrixde,
    by.x = "Gene.stable.ID",
    by.y = "row.names",
    all.X = T, all.Y = T
)

ExpressMatrixde_HGNC <- ExpressMatrixde_HGNC[, !(names(ExpressMatrixde_HGNC) %in% c("Gene.stable.ID"))]
ExpressMatrixde_HGNC <- avereps(ExpressMatrixde_HGNC,
    ID = ExpressMatrixde_HGNC$HGNC.symbol
)
rownames(ExpressMatrixde_HGNC) <- ExpressMatrixde_HGNC[, "HGNC.symbol"]
ExpressMatrixde_HGNC <- ExpressMatrixde_HGNC[, !(colnames(ExpressMatrixde_HGNC) %in% c("HGNC.symbol"))]
ExpressMatrixde_HGNC <- as.matrix(data.frame(ExpressMatrixde_HGNC, check.names = FALSE))
class(ExpressMatrixde_HGNC) <- "numeric"
write.csv(ExpressMatrixde_HGNC, file.path(deresults_path, "expression_matrix_de_lm_hgnc.csv"), quote = F)

dataMatrixde_HGNC <- merge(rhesus2human, dataMatrixde,
    by.x = "Gene.stable.ID",
    by.y = "row.names",
    all.X = T, all.Y = T
)

dataMatrixde_HGNC <- dataMatrixde_HGNC[, !(names(dataMatrixde_HGNC) %in% c("Gene.stable.ID"))]
dataMatrixde_HGNC <- avereps(dataMatrixde_HGNC,
    ID = dataMatrixde_HGNC$HGNC.symbol
)
rownames(dataMatrixde_HGNC) <- dataMatrixde_HGNC[, "HGNC.symbol"]
dataMatrixde_HGNC <- dataMatrixde_HGNC[, !(colnames(dataMatrixde_HGNC) %in% c("HGNC.symbol"))]
dataMatrixde_HGNC <- as.matrix(data.frame(dataMatrixde_HGNC, check.names = FALSE))
class(dataMatrixde_HGNC) <- "numeric"
write.csv(dataMatrixde_HGNC, file.path(deresults_path, "full_expression_matrix_de_hgnc.csv"), quote = F)

message(paste0("Dimensionality of DE genes ", dim(ExpressMatrixde)[1]))
new_colnames <- rep(c("D3", "D14", "D105", "D108", "D119"), 4)
colnames(ExpressMatrixde) <- new_colnames

new_colnames <- c(
    c("68.1.CyCMV_D3.Prot", "68.1.CyCMV_D14.Prot", "68.1.CyCMV_D105.Prot", "68.1.CyCMV_D108.Prot", "68.1.CyCMV_D119.Prot"),
    c("68.1.CyCMV_D3.NonProt", "68.1.CyCMV_D14.NonProt", "68.1.CyCMV_D105.NonProt", "68.1.CyCMV_D108.NonProt", "68.1.CyCMV_D119.NonProt"),
    c("full.length.CyCMV_D3.Prot", "full.length.CyCMV_D14.Prot", "full.length.CyCMV_D105.Prot", "full.length.CyCMV_D108.Prot", "full.length.CyCMV_D119.Prot"),
    c("full.length.CyCMV_D3.NonProt", "full.length.CyCMV_D14.NonProt", "full.length.CyCMV_D105.NonProt", "full.length.CyCMV_D108.NonProt", "full.length.CyCMV_D119.NonProt")
)

colnames(results) <- new_colnames
colnames(dataMatrixde) <- new_colnames
colnames(dataMatrixde_HGNC) <- new_colnames
colnames(ExpressMatrixde_HGNC) <- new_colnames

# --barplot
vizualize_DE_genes_bp(results, file.path(deresults_path, "barplot.svg"))

# --heatmap
colcolorlistprot <- rep(c(rep("red", 5), rep("black", 5)), 2)
colcolorvaccine <- c(rep("#068d06", 10), rep("orange", 10))
colcolormatrix <- cbind(colcolorlistprot, colcolorvaccine)
colcolormatrix <- as.matrix(colcolormatrix)
colnames(colcolormatrix) <- c("Protection", "Vaccine")

png(file.path(deresults_path, "heatmap.png"), width = 7, height = 8, units = "in", res=300)
# svglite(file.path(deresults_path, "heatmap.svg"), width = 7, height = 8)
global_modulesde <- heatmap.L.4(ExpressMatrixde,
    figmargins = c(7, 5), colcolorlist=colcolormatrix,
    cutoff = 1, distmethod = "pearson", cexcol = 1, colsep=10, sepcolor="white",
    clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 1.4, sepwidth=0.2
)
dev.off()


######### -- Signature Overlap --#########
results_folder <- "Figure3_heatmap"
generate_folder(results_folder)

#### -- generate heatmap of Dan's 122 genes from manuscript with Just Cyno LFCs -- ####
dangenes <- read.csv("./Maloulietal2022_SigGenes.csv", row.names = 1)
dde <- read.csv("Barranesetal2021_SigGenes.csv", header = TRUE, row.names = 1, quote = "")

dangenes <- dangenes[, 1:18]
genes <- intersect(rownames(dde), rownames(dataMatrixde))

cymvdde <- dataMatrixde_HGNC[rownames(dataMatrixde_HGNC) %in% rownames(dangenes), ]
cymvdde <- merge(dangenes, cymvdde, by.x = "row.names", by.y = "row.names")
rownames(cymvdde) <- cymvdde$Row.names
cymvdde$Row.names <- NULL
message("STATUS dimensions of all Maloulietal2022 genes ", dim(cymvdde)[1])


colcolorlistprot <- c(rep(c(rep("red", 9), rep("black", 9)), 1), rep(c(rep("red", 5), rep("black", 5)), 2))
colcolorvaccine <- c(rep("gray", 18), rep("#068d06", 10), rep("orange", 10))
colcolormatrix <- cbind(colcolorlistprot, colcolorvaccine)
colcolormatrix <- as.matrix(colcolormatrix)
colnames(colcolormatrix) <- c("Protection", "Vaccine")

cymvdde <- as.matrix(as.data.frame(cymvdde))
class(cymvdde) <- "numeric"

png(file.path(results_folder, "heatmap.png"), width = 9, height = 8, units = "in", res = 300)
# svglite(file.path(results_folder, "heatmap_alldde.svg"), width = 9, height = 8)
global_modulesdeovall <- heatmap.L.4(cymvdde,
    figmargins = c(7, 5), colcolorlist = colcolormatrix,
    cutoff = 1, distmethod = "pearson", cexcol = 1, 
    colsep = c(9, 18, 23, 28, 33, 38),
    sepcolor = "white", sepwidth = 0.05,
    clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 1.4,
    labCol = FALSE
)
dev.off()
totaldata <- merge(as.data.frame(global_modulesdeovall$modulesrows), 
    global_modulesdeovall$clustermatrix, by.x="row.names", by.y="row.names")
rownames(totaldata)<-totaldata$Row.names
totaldata$Row.names <- NULL
colnames(totaldata)[1] <- "cluster"
totaldata <- totaldata[rownames(global_modulesdeovall$clustermatrix),]
write.csv(totaldata, file.path(results_folder, "totaldata4hm.csv"))

