# Calculate RareEffect PRS for each individual

## Load packages and sources

library(SAIGE)
library(data.table)
library(Matrix)
library(optparse)

option_list <- list(
    make_option(c("--effectFile"), type="character", default="",
        help="Path to effect file from RareEffect main function"),
    make_option(c("--bedFile"), type="character", default="",
        help="Path to bed file containing rare variants"),
    make_option(c("--bimFile"), type="character", default="",
        help="Path to bim file containing rare variants"),
    make_option(c("--famFile"), type="character", default="",
        help="Path to fam file containing rare variants"),
    make_option(c("--groupFile"), type="character", default="",
        help="Path to group file (containing functional annotation of variants)"),
    make_option(c("--geneName"), type="character", default="",
        help="Gene name to analyze"),
    make_option(c("--commonVariantsToExclude"), type="character", default="",
        help="Path to file containing list of common variant IDs to exclude from PRS calculation (one variant ID per line)"),
    make_option(c("--outputPrefix"), type="character", default="",
        help="Path to save output (without extension and gene name)")
)

## list of options
parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = FALSE)
print(args)

## Validate required inputs
required_args <- c("effectFile", "bedFile", "bimFile", "famFile", "groupFile", "geneName", "outputPrefix")
for (arg_name in required_args) {
    if (args[[arg_name]] == "") {
        stop(paste0("--", arg_name, " is required but not provided."))
    }
}
if (!file.exists(args$effectFile)) stop(paste0("Effect file not found: ", args$effectFile))
if (!file.exists(args$bedFile)) stop(paste0("Bed file not found: ", args$bedFile))

## Assign variables
bedFile <- args$bedFile
bimFile <- args$bimFile
famFile <- args$famFile
effectFile <- args$effectFile
groupFile <- args$groupFile
geneName <- args$geneName

## Set PLINK object
bim <- fread(bimFile)
fam <- fread(famFile)
sampleID <- as.character(fam$V2)
modglmm <- list()
modglmm$sampleID <- sampleID
n_samples <- length(sampleID)
var_by_func_anno <- read_groupfile(groupFile, geneName)

objGeno <- SAIGE::setGenoInput(bgenFile = "",
        bgenFileIndex = "",
        vcfFile = "",
        vcfFileIndex = "",
        vcfField = "",
        savFile = "",
        savFileIndex = "",
        sampleFile = "",
        bedFile=bedFile,
        bimFile=bimFile,
        famFile=famFile,
        idstoIncludeFile = "",
        rangestoIncludeFile = "",
        chrom = "",
        AlleleOrder = "alt-first",
        sampleInModel = sampleID
    )

## Read effect file and restore minor-allele-relative effects for flipped variants
effect <- fread(effectFile)
if ("is_flipped" %in% colnames(effect)) {
    effect[is_flipped == TRUE, effect := -effect]
}

# LoF
if (length(var_by_func_anno[[1]]) == 0) {
    print("LoF variant does not exist.")
}
# missense
if (length(var_by_func_anno[[2]]) == 0) {
    print("missense variant does not exist.")
}
# synonymous
if (length(var_by_func_anno[[3]]) == 0) {
    print("synonymous variant does not exist.")
}

# Remove variants not in plink file
for (i in 1:3) {
    var_by_func_anno[[i]] <- var_by_func_anno[[i]][which(var_by_func_anno[[i]] %in% objGeno$markerInfo$ID)]
}
print(str(var_by_func_anno))

# Read genotype matrix
lof_mat_collapsed_all <- collapse_matrix(objGeno, var_by_func_anno[[1]], sampleID, modglmm, macThreshold = 0)
lof_mat_collapsed <- lof_mat_collapsed_all[[1]]
lof_flipped <- lof_mat_collapsed_all[[2]]

mis_mat_collapsed_all <- collapse_matrix(objGeno, var_by_func_anno[[2]], sampleID, modglmm, macThreshold = 0)
mis_mat_collapsed <- mis_mat_collapsed_all[[1]]
mis_flipped <- mis_mat_collapsed_all[[2]]

syn_mat_collapsed_all <- collapse_matrix(objGeno, var_by_func_anno[[3]], sampleID, modglmm, macThreshold = 0)
syn_mat_collapsed <- syn_mat_collapsed_all[[1]]
syn_flipped <- syn_mat_collapsed_all[[2]]
print("Genotype matrix is loaded.")

# Remove common variants
if (args$commonVariantsToExclude == "") {
    print("No common variant exclusion file is provided.")
} else {
    var_list <- fread(args$commonVariantsToExclude, header = FALSE)
    common_var <- var_list$V1
    lof_mat_collapsed <- lof_mat_collapsed[, !colnames(lof_mat_collapsed) %in% common_var, drop = FALSE]
    mis_mat_collapsed <- mis_mat_collapsed[, !colnames(mis_mat_collapsed) %in% common_var, drop = FALSE]
    syn_mat_collapsed <- syn_mat_collapsed[, !colnames(syn_mat_collapsed) %in% common_var, drop = FALSE]
    print(paste0("Common variants removed. Remaining columns - lof: ", ncol(lof_mat_collapsed),
                 ", mis: ", ncol(mis_mat_collapsed), ", syn: ", ncol(syn_mat_collapsed)))
}

## Helper: map effect sizes from effect file to genotype matrix columns
map_effects <- function(mat, effect_df, ur_label) {
    nc <- ncol(mat)
    if (nc == 0) {
        return(numeric(0))
    }

    # Get UR effect (default 0 if not in effect file)
    ur_idx <- which(effect_df$variant == ur_label)
    ur_effect <- if (length(ur_idx) > 0) effect_df$effect[ur_idx] else 0

    # Initialize all positions with UR effect
    eff_vec <- rep(ur_effect, nc)

    # Overwrite with individual variant effects where available
    cnames <- colnames(mat)
    for (i in seq_len(nc)) {
        idx <- which(effect_df$variant == cnames[i])
        if (length(idx) > 0) {
            eff_vec[i] <- effect_df$effect[idx]
        }
    }
    return(eff_vec)
}

# Read effect size
effect <- as.data.frame(effect)

effect_lof <- map_effects(lof_mat_collapsed, effect, "lof_UR")
effect_mis <- map_effects(mis_mat_collapsed, effect, "mis_UR")
effect_syn <- map_effects(syn_mat_collapsed, effect, "syn_UR")

print("Effect size is loaded.")

## Calculate PRS
calc_group_prs <- function(mat, eff) {
    if (ncol(mat) == 0) return(Matrix::Matrix(0, nrow = nrow(mat), ncol = 1, sparse = TRUE))
    mat %*% eff
}

lof_prs <- calc_group_prs(lof_mat_collapsed, effect_lof)
mis_prs <- calc_group_prs(mis_mat_collapsed, effect_mis)
syn_prs <- calc_group_prs(syn_mat_collapsed, effect_syn)
total_prs <- lof_prs + mis_prs + syn_prs

out <- data.frame(
    IID = rownames(total_prs),
    PRS = as.numeric(total_prs[, 1]),
    stringsAsFactors = FALSE
)

## Save PRS
write.table(out, paste0(args$outputPrefix, "_", geneName, "_score.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
