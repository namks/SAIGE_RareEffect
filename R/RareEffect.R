# Read Group file and split variants by functional annotations
read_groupfile <- function(groupfile_name, gene_name) {
    lines <- readLines(groupfile_name)

    # Pre-filter to only the gene's lines via fixed-string startsWith
    # (regex strsplit on every line is the dominant cost for large group files)
    prefix_sp <- paste0(gene_name, " ")
    prefix_tab <- paste0(gene_name, "\t")
    gene_lines <- lines[startsWith(lines, prefix_sp) | startsWith(lines, prefix_tab)]

    var <- NULL
    anno <- NULL
    for (marker_group_line in gene_lines) {
        marker_group_line_list <- strsplit(marker_group_line, split = "[ \t]+")[[1]]
        if (length(marker_group_line_list) < 2) next
        if (marker_group_line_list[2] == "var") {
            var <- marker_group_line_list
        } else {
            anno <- marker_group_line_list
        }
        if (!is.null(var) && !is.null(anno)) break
    }

    lof_var <- var[which(anno == "lof")]
    mis_var <- var[which(anno == "missense")]
    syn_var <- var[which(anno == "synonymous")]
    list(lof_var, mis_var, syn_var)
}

read_matrix_by_one_marker <- function(objGeno, var_list, sampleID) {
    n_samples <- length(sampleID)
    if (length(var_list) == 0) {
        mat <- Matrix::Matrix(0, nrow = n_samples, ncol = 0, sparse = TRUE)
        message("No variants in the list")
        return(list(mat, rep(FALSE, 0)))
    }
    p <- length(var_list)
    is_flipped <- rep(FALSE, p)
    i_list <- vector("list", p)
    x_list <- vector("list", p)
    nz_counts <- integer(p)

    # Precompute marker index in objGeno$markerInfo$ID once
    idx_map <- match(var_list, objGeno$markerInfo$ID)
    if (anyNA(idx_map)) {
        missing_vars <- var_list[is.na(idx_map)]
        stop(sprintf("Variants not found in objGeno$markerInfo$ID: %s%s",
                     paste(head(missing_vars, 5), collapse = ", "),
                     if (length(missing_vars) > 5) sprintf(" (+%d more)", length(missing_vars) - 5) else ""))
    }

    for (i in seq_along(var_list)) {
        t_GVec <- rep(0, n_samples)
        idx <- idx_map[i]
        SAIGE::Unified_getOneMarker(t_genoType = objGeno$genoType,
            t_gIndex_prev = objGeno$markerInfo$genoIndex_prev[idx],
            t_gIndex = objGeno$markerInfo$genoIndex[idx],
            t_ref = "2",
            t_alt = "1",
            t_marker = objGeno$markerInfo$ID[idx],
            t_pd = objGeno$markerInfo$POS[idx],
            t_chr = toString(objGeno$markerInfo$CHROM[idx]),
            t_altFreq = 0,
            t_altCounts = 0,
            t_missingRate = 0,
            t_imputeInfo = 0,
            t_isOutputIndexForMissing = TRUE,
            t_indexForMissing = 0,
            t_isOnlyOutputNonZero = FALSE,
            t_indexForNonZero = 0,
            t_GVec = t_GVec,
            t_isImputation = FALSE
        )

        # Single scan: positive entries are the non-missing, non-zero genotypes
        # (missing is encoded as -1 by C++ side and is skipped naturally)
        nz <- which(t_GVec > 0)
        vals <- t_GVec[nz]
        s <- sum(vals)
        if (s > n_samples) {
            # Major allele is ALT — flip. Almost never hit for rare variants.
            # Note: this densifies the column (most zeros become 2s).
            t_GVec[t_GVec < 0] <- 0
            t_GVec <- 2 - t_GVec
            is_flipped[i] <- TRUE
            nz <- which(t_GVec > 0)
            vals <- t_GVec[nz]
        }
        i_list[[i]] <- nz
        x_list[[i]] <- vals
        nz_counts[i] <- length(nz)
    }

    # Construct sparseMatrix from accumulated triplets in one call
    mat <- Matrix::sparseMatrix(
        i = unlist(i_list, use.names = FALSE),
        j = rep.int(seq_len(p), nz_counts),
        x = unlist(x_list, use.names = FALSE),
        dims = c(n_samples, p),
        dimnames = list(sampleID, var_list)
    )

    return(list(mat, is_flipped))
}

collapse_matrix <- function(objGeno, var_list, sampleID, modglmm, macThreshold = 10) {
    n_samples <- length(sampleID)
    mat <- Matrix::Matrix(0, nrow = n_samples, ncol = 0, sparse = TRUE)
    if (length(var_list) == 0) {
        message("No variants in the list")
        return(list(mat, NULL))
    }
    tmp <- read_matrix_by_one_marker(objGeno, var_list, sampleID)
    mat <- tmp[[1]]
    is_flipped <- tmp[[2]]
    flipped_var <- as.data.table(cbind(var_list, is_flipped))
    mat <- mat[which(rownames(mat) %in% modglmm$sampleID), , drop = FALSE]
    MAC <- colSums(mat)
    MAF <- MAC / (2 * nrow(mat))
    idx_rare <- which((MAF < 0.01) & (MAC >= macThreshold))
    idx_UR <- which((MAF > 0) & (MAC < macThreshold))
    if (length(idx_rare) == 0 & length(idx_UR) == 0) {
        message("No variants with MAF < 0.01 or MAC > 0")
        empty_mat <- Matrix::Matrix(0, nrow = nrow(mat), ncol = 0, sparse = TRUE)
        return(list(empty_mat, flipped_var))
    }
    mat_rare <- mat[, idx_rare, drop = FALSE]
    mat_UR <- mat[, idx_UR, drop = FALSE]
    UR_rowsum <- rowSums(mat_UR)
    UR_rowsum[which(UR_rowsum > 1)] <- 1
    mat_UR_collapsed <- cbind(mat_rare, UR_rowsum)
    return(list(mat_UR_collapsed, flipped_var))
}

get_range <- function(v) {
    # input are vectors of variants by functional annotation
    start_pos <- Inf
    end_pos <- 0
    for (i in 1:3) {
        if (length(v[[i]]) > 0) {
            pos <- as.numeric(str_split_fixed(v[[i]], ":", 4)[,2]) # Modified in 0.3
            sub_start_pos <- min(pos)
            sub_end_pos <- max(pos)
            if (start_pos > sub_start_pos) {
                start_pos <- sub_start_pos
            }
            if (end_pos < sub_end_pos) {
                end_pos <- sub_end_pos
            }
        }
    }
    return (list(start_pos, end_pos))
}


# Read phenotype file
read_pheno <- function(pheno_file, pheno_code, iid_col = "f.eid") {
    pheno <- data.table::fread(pheno_file, quote = "")
    out <- subset(pheno, select = c(iid_col, pheno_code))
    return(out)
}


calc_log_lik <- function(delta, S, UtY, Y_UUtY) {
    k <- length(S)
    n <- nrow(Y_UUtY)

    log_lik1 <- 0
    for (i in 1:k) {
        log_lik1 <- log_lik1 + (UtY[i, ])^2 / (S[i] + delta)
    }

    log_lik2 <- 1 / delta * sum((Y_UUtY)^2)

    out <- -0.5 * (n * log(2 * pi) + sum(log(S + delta)) + (n - k) * log(delta)
                   + n + n * log(1 / n * (log_lik1 + log_lik2)))

    return(as.numeric(out))
}

calc_post_beta <- function(G, delta, S, UtY, U, Sigma = NULL) {
    d_vec <- as.numeric(1 / (S + delta))
    if (length(d_vec) == 1) {
        D <- matrix(d_vec, 1, 1)
    } else {
        D <- diag(d_vec)
    }
    UDUtY <- U %*% D %*% UtY
    if (is.null(Sigma)) {
        out <- t(G) %*% UDUtY
    } else {
        K_sparse <- as(Sigma, "dgCMatrix")
        out <- K_sparse %*% t(G) %*% UDUtY
    }
    return(out)
}

# Run FaST-LMM to obtain posterior beta
fast_lmm <- function(G, Y, Sigma = NULL) {
    G[is.na(G)] <- 0
    GtG <- crossprod(G)
    if (sum(G, na.rm = TRUE) == 0) {
        return(list(as.matrix(0), 0, 1e6, GtG))
    }
    Y <- as.matrix(Y)

    # W = G %*% chol(Sigma), or just G when Sigma = I
    if (is.null(Sigma)) {
        W_sparse <- as(G, "dgCMatrix")
    } else {
        L <- chol(Sigma)
        W_sparse <- as(G %*% L, "dgCMatrix")
    }

    svd_mat <- sparsesvd::sparsesvd(W_sparse)
    U <- svd_mat$u
    S <- (svd_mat$d)^2
    UtY <- t(U) %*% Y
    Y_UUtY <- Y - U %*% UtY

    opt <- optim(par = 1, fn = calc_log_lik, S = S, UtY = UtY, Y_UUtY = Y_UUtY,
                 method = "Brent", lower = 1e-10, upper = 1e6,
                 control = list(fnscale = -1))

    tr_GtG_Sigma <- if (is.null(Sigma)) sum(diag(GtG)) else sum(diag(GtG %*% Sigma))
    post_beta <- calc_post_beta(G, opt$par, S, UtY, U, Sigma = Sigma)

    return(list(post_beta, tr_GtG_Sigma, opt$par, GtG))
}


calc_gene_effect_size <- function(G, lof_ncol, post_beta, Y) {
    beta_lof <- post_beta[1:lof_ncol]
    beta_lof <- abs(beta_lof)
    lof_prs <- G[, 1:lof_ncol, drop = F] %*% beta_lof
    lof_prs_norm <- (lof_prs - mean(lof_prs)) / sd(lof_prs)
    m1 <- lm(Y ~ lof_prs_norm[, 1])

    return(m1$coefficients[2])
}

mom_estimator_marginal <- function(G, y, Sigma = NULL) {
    n <- length(y)
    G <- as(G, "dgCMatrix")
    G[is.na(G)] <- 0

    if (is.null(Sigma)) {
        GtG <- crossprod(G)
        t1 <- sum(GtG^2)                    # tr((GG')^2)
        t2 <- sum(G^2)                      # tr(G'G)
        Gty <- crossprod(G, y)
        c1 <- as.numeric(sum(Gty^2))        # y'GG'y
    } else {
        GtG_Sigma <- crossprod(G) %*% Sigma
        t1 <- sum(GtG_Sigma^2)
        t2 <- sum(diag(GtG_Sigma))
        c1 <- as.numeric(t(y) %*% G %*% Sigma %*% t(G) %*% y)
    }

    c2 <- sum(y^2)
    A <- matrix(c(t1, t2, t2, n), ncol = 2)
    b <- matrix(c(c1, c2), ncol = 1)
    return(solve(A, b))
}

mom_estimator_joint <- function(G1, G2, G3, y,
                                 Sigma1 = NULL, Sigma2 = NULL, Sigma3 = NULL) {
    n <- length(y)
    G1 <- as(G1, "dgCMatrix"); G1[is.na(G1)] <- 0
    G2 <- as(G2, "dgCMatrix"); G2[is.na(G2)] <- 0
    G3 <- as(G3, "dgCMatrix"); G3[is.na(G3)] <- 0

    # W_k = G_k %*% chol(Sigma_k), or G_k when Sigma_k = I
    W1 <- if (is.null(Sigma1)) G1 else G1 %*% chol(Sigma1)
    W2 <- if (is.null(Sigma2)) G2 else G2 %*% chol(Sigma2)
    W3 <- if (is.null(Sigma3)) G3 else G3 %*% chol(Sigma3)

    t11 <- sum(crossprod(W1)^2)
    t22 <- sum(crossprod(W2)^2)
    t33 <- sum(crossprod(W3)^2)
    t12 <- sum(crossprod(W1, W2)^2)
    t13 <- sum(crossprod(W1, W3)^2)
    t23 <- sum(crossprod(W2, W3)^2)
    t14 <- sum(W1^2)
    t24 <- sum(W2^2)
    t34 <- sum(W3^2)

    A <- matrix(c(t11,t12,t13,t14, t12,t22,t23,t24,
                  t13,t23,t33,t34, t14,t24,t34,n), ncol = 4)

    c1 <- as.numeric(sum(crossprod(W1, y)^2))
    c2 <- as.numeric(sum(crossprod(W2, y)^2))
    c3 <- as.numeric(sum(crossprod(W3, y)^2))
    b <- matrix(c(c1, c2, c3, sum(y^2)), ncol = 1)

    return(solve(A, b))
}

# calculate_single_blup <- function(G, delta, Sigma, y) {
#     # Henderson Mixed Model Equation: beta = (G^T G + delta * Sigma^-1)^-1 * G^T y
#     G <- as(G, "dgCMatrix")
#     beta <- solve(t(G) %*% G + delta * solve(Sigma)) %*% t(G) %*% y
#     return (beta)
# }

calculate_joint_blup <- function(G1, G2, G3, tau1, tau2, tau3, psi, y) {
    G1 <- as(G1, "dgCMatrix")
    G2 <- as(G2, "dgCMatrix")
    G3 <- as(G3, "dgCMatrix")
    G <- cbind(G1, G2, G3)
    G[is.na(G)] <- 0
    p1 <- ncol(G1); p2 <- ncol(G2); p3 <- ncol(G3)

    GtG <- crossprod(G)
    Gty <- crossprod(G, y)

    # (G'G / psi + Sigma^{-1}) where Sigma = bdiag(tau1*I, tau2*I, tau3*I)
    C <- GtG / psi
    inv_tau_diag <- c(rep(1/tau1, p1), rep(1/tau2, p2), rep(1/tau3, p3))
    diag(C) <- diag(C) + inv_tau_diag

    beta <- solve(C, Gty / psi)
    return(beta)
}

weight_cal <- function(beta_k, delta = 10^(-5), gamma = 2, q = 0, factor2 = 1) {
	p <- length(beta_k)
	w_out <- rep(0, p)
	idx <- which(abs(beta_k) <= delta)
	if(length(idx)> 0){
		beta_k1<-beta_k[idx]
		b_delta = abs(beta_k1/delta)^{gamma}
		logw1 = (q-2) * log(delta) + (q-2)/gamma* log(1+ b_delta)
		w_out[idx]<-exp(logw1 * factor2)
	}
	idx<-which(abs(beta_k) > delta)
	if(length(idx)> 0){
		beta_k1<-beta_k[idx]
		b_delta_inv = abs(delta/beta_k1)^{gamma}
		logw2 = (q-2) * log(abs(beta_k1)) + (q-2)/gamma* log(1+ b_delta_inv)
		w_out[idx]<-exp(logw2)
	}

	# make sum of all to be p
	a1 = 1/w_out
	a1_out = a1/sum(a1)*p
	w_out = 1/a1_out

	return(list(w_out=w_out))
}

adaptive_ridge <- function(X, y, lambda, q = 0, delta = 1e-5, gamma = 2, max_iter = 100, tol = 0.01, sigma_sq = 1) {
    w <- rep(1, ncol(X))
    beta <- rep(0, ncol(X))

    Xy <- crossprod(X, y)
    XX <- crossprod(X)
    for (k in 1:max_iter) {
        A <- XX + 0                                 # copy to preserve XX
        diag(A) <- diag(A) + lambda * w
        beta_new <- solve(A, Xy)

        w_new <- weight_cal(beta_new, delta = delta, gamma = gamma, q = q)$w_out

        if (sqrt(sum((beta_new - beta)^2)) < sum(beta^2) * tol) {
            break
        }

        beta <- beta_new
        w <- w_new
    }

    A <- XX + 0
    diag(A) <- diag(A) + lambda * w
    A_inv <- chol2inv(chol(A))

    cov_beta <- as.numeric(sigma_sq) * (A_inv %*% XX %*% A_inv)
    se_beta <- sqrt(diag(cov_beta))

    list(beta = beta, weights = w, iterations = k, se_beta = se_beta)
}

# Helper: run FaST-LMM + MoM for a single annotation group
estimate_group <- function(mat, y_vec, sigma_sq) {
    result <- fast_lmm(G = mat, Y = y_vec)
    post_beta <- result[[1]]
    tr_GtG    <- result[[2]]
    delta     <- result[[3]]
    GtG       <- result[[4]]
    tau       <- as.numeric(sigma_sq / delta)
    tau_mom   <- mom_estimator_marginal(G = mat, y = y_vec)
    list(post_beta = post_beta, tr_GtG = tr_GtG, delta = delta,
         GtG = GtG, tau = tau, tau_mom = tau_mom)
}

# Helper: compute gene-level heritability for one group
compute_h2 <- function(tau_adj, tr_GtG, sigma_sq, denom_extra) {
    max(tau_adj * tr_GtG / (tau_adj * tr_GtG + sigma_sq * denom_extra), 0)
}

# Helper: compute PEV for one group
compute_pev <- function(GtG, sigma_sq, tau_adj) {
    GtG_scaled <- GtG / as.numeric(sigma_sq)
    diag(GtG_scaled) <- diag(GtG_scaled) + 1 / tau_adj
    diag(solve(GtG_scaled))
}

run_RareEffect <- function(rdaFile, chrom, geneName, groupFile, traitType, bedFile, bimFile, famFile, macThreshold, collapseLoF, collapsemis, collapsesyn, apply_AR, outputPrefix) {
    # Load SAIGE step 1 results
    load(rdaFile)

    if (traitType == "binary") {
        # modglmm$residuals <- modglmm$Y - modglmm$linear.predictors
        v <- modglmm$fitted.values * (1 - modglmm$fitted.values)    # v_i = mu_i * (1 - mu_i)
        modglmm$residuals <- sqrt(v) * (1 / (modglmm$fitted.values * (1 - modglmm$fitted.values)) * (modglmm$y - modglmm$fitted.values))
    }
    sigma_sq <- var(modglmm$residuals)

    # Set PLINK object
    # Only column 2 of fam (IID) is needed; bim is not used at R level
    fam <- fread(famFile, select = 2, header = FALSE)
    sampleID <- as.character(fam[[1]])
    n_samples <- length(sampleID)

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

    message("Analysis started")

    var_by_func_anno <- read_groupfile(groupFile, geneName)

    # LoF
    if (length(var_by_func_anno[[1]]) == 0) {
        message("LoF variant does not exist.")
    }
    # missense
    if (length(var_by_func_anno[[2]]) == 0) {
        message("missense variant does not exist.")
    }
    # synonymous
    if (length(var_by_func_anno[[3]]) == 0) {
        message("synonymous variant does not exist.")
    }

    # Remove variants not in plink file
    for (i in 1:3) {
        var_by_func_anno[[i]] <- var_by_func_anno[[i]][which(var_by_func_anno[[i]] %in% objGeno$markerInfo$ID)]
    }
    # print(str(var_by_func_anno))

    # Read genotype matrix
    if (collapseLoF) {
        lof_mat_collapsed_all <- collapse_matrix(objGeno, var_by_func_anno[[1]], sampleID, modglmm, macThreshold = 0)
        lof_mat_collapsed <- lof_mat_collapsed_all[[1]]
        lof_flipped <- lof_mat_collapsed_all[[2]]
        lof_mat_collapsed <- Matrix::Matrix(rowSums(lof_mat_collapsed), ncol = 1, sparse = TRUE, dimnames = list(rownames(lof_mat_collapsed), NULL))
    } else {
        lof_mat_collapsed_all <- collapse_matrix(objGeno, var_by_func_anno[[1]], sampleID, modglmm, macThreshold)
        lof_mat_collapsed <- lof_mat_collapsed_all[[1]]
        lof_flipped <- lof_mat_collapsed_all[[2]]
    }
    if (collapsemis) {
        mis_mat_collapsed_all <- collapse_matrix(objGeno, var_by_func_anno[[2]], sampleID, modglmm, macThreshold = 0)
        mis_mat_collapsed <- mis_mat_collapsed_all[[1]]
        mis_flipped <- mis_mat_collapsed_all[[2]]
        mis_mat_collapsed <- Matrix::Matrix(rowSums(mis_mat_collapsed), ncol = 1, sparse = TRUE, dimnames = list(rownames(mis_mat_collapsed), NULL))
    } else {
        mis_mat_collapsed_all <- collapse_matrix(objGeno, var_by_func_anno[[2]], sampleID, modglmm, macThreshold)
        mis_mat_collapsed <- mis_mat_collapsed_all[[1]]
        mis_flipped <- mis_mat_collapsed_all[[2]]
    }

    if (collapsesyn) {
        syn_mat_collapsed_all <- collapse_matrix(objGeno, var_by_func_anno[[3]], sampleID, modglmm, macThreshold = 0)
        syn_mat_collapsed <- syn_mat_collapsed_all[[1]]
        syn_flipped <- syn_mat_collapsed_all[[2]]
        syn_mat_collapsed <- Matrix::Matrix(rowSums(syn_mat_collapsed), ncol = 1, sparse = TRUE, dimnames = list(rownames(syn_mat_collapsed), NULL))
    } else {
        syn_mat_collapsed_all <- collapse_matrix(objGeno, var_by_func_anno[[3]], sampleID, modglmm, macThreshold)
        syn_mat_collapsed <- syn_mat_collapsed_all[[1]]
        syn_flipped <- syn_mat_collapsed_all[[2]]
    }

    # Set column name
    if (ncol(lof_mat_collapsed) > 0) {
        colnames(lof_mat_collapsed)[ncol(lof_mat_collapsed)] <- "lof_UR"
    }
    if (ncol(mis_mat_collapsed) > 0) {
        colnames(mis_mat_collapsed)[ncol(mis_mat_collapsed)] <- "mis_UR"
    }
    if (ncol(syn_mat_collapsed) > 0) {
        colnames(syn_mat_collapsed)[ncol(syn_mat_collapsed)] <- "syn_UR"
    }

    # Make genotype matrix
    nonempty_mat <- list()
    if (ncol(lof_mat_collapsed) > 0) nonempty_mat <- c(nonempty_mat, list(lof_mat_collapsed))
    if (ncol(mis_mat_collapsed) > 0) nonempty_mat <- c(nonempty_mat, list(mis_mat_collapsed))
    if (ncol(syn_mat_collapsed) > 0) nonempty_mat <- c(nonempty_mat, list(syn_mat_collapsed))

    G <- do.call(cbind, nonempty_mat)
    lof_ncol <- ncol(lof_mat_collapsed)
    mis_ncol <- ncol(mis_mat_collapsed)
    syn_ncol <- ncol(syn_mat_collapsed)

    # print("Dimension of G")
    # print(dim(G))

    # Obtain residual vector and genotype matrix with the same sample order
    # Use index-based alignment instead of cbind (which coerces numeric to character)
    common_samples <- intersect(modglmm$sampleID, rownames(G))
    model_idx <- match(common_samples, modglmm$sampleID)
    y_vec <- as.numeric(modglmm$residuals[model_idx])
    G_reordered <- G[match(common_samples, rownames(G)), , drop = FALSE]
    if (traitType == "binary") {
        v_ordered <- v[model_idx]
        vG_reordered <- as.vector(sqrt(v_ordered)) * G_reordered    # Sigma_e^(-1/2) G
    }
    n_samples <- nrow(G_reordered)

    # Define matrices by functional annotation
    if (lof_ncol == 0) {
        lof_mat <- NULL
    } else {
        if (traitType == "binary") {
            lof_mat <- vG_reordered[,c(1:lof_ncol), drop = F]
        } else {
            lof_mat <- G_reordered[,c(1:lof_ncol), drop = F]
        }
    }

    if (mis_ncol == 0) {
        mis_mat <- NULL
    } else {
        if (traitType == "binary") {
            mis_mat <- vG_reordered[,c((lof_ncol + 1):(lof_ncol + mis_ncol)), drop = F]
        } else {
            mis_mat <- G_reordered[,c((lof_ncol + 1):(lof_ncol + mis_ncol)), drop = F]
        }
    }

    if (syn_ncol == 0) {
        syn_mat <- NULL
    } else {
        if (traitType == "binary") {
            syn_mat <- vG_reordered[,c((lof_ncol + mis_ncol + 1):(lof_ncol + mis_ncol + syn_ncol)), drop = F]
        } else {
            syn_mat <- G_reordered[,c((lof_ncol + mis_ncol + 1):(lof_ncol + mis_ncol + syn_ncol)), drop = F]
        }
    }

    # Run FaST-LMM + MoM for each functional annotation group
    est_lof <- if (lof_ncol > 0) estimate_group(lof_mat, y_vec, sigma_sq) else NULL
    est_mis <- if (mis_ncol > 0) estimate_group(mis_mat, y_vec, sigma_sq) else NULL
    est_syn <- if (syn_ncol > 0) estimate_group(syn_mat, y_vec, sigma_sq) else NULL

    tau_lof <- if (!is.null(est_lof)) est_lof$tau else 0
    tau_mis <- if (!is.null(est_mis)) est_mis$tau else 0
    tau_syn <- if (!is.null(est_syn)) est_syn$tau else 0

    # Estimate variance component tau jointly (MoM estimator), only if all groups are non-empty
    if ((lof_ncol > 0) & (mis_ncol > 0) & (syn_ncol > 0)) {
        tau_mom_joint <- mom_estimator_joint(
            G1 = lof_mat,
            G2 = mis_mat,
            G3 = syn_mat,
            y = y_vec
        )

        # Adjust variance component tau
        # Only apply ratio correction when both marginal and joint MoM estimates are positive
        if (est_lof$tau_mom[1] > 0 && tau_mom_joint[1] > 0) {
            tau_lof_adj <- tau_lof * tau_mom_joint[1] / est_lof$tau_mom[1]
        } else {
            tau_lof_adj <- tau_lof
        }
        if (est_mis$tau_mom[1] > 0 && tau_mom_joint[2] > 0) {
            tau_mis_adj <- tau_mis * tau_mom_joint[2] / est_mis$tau_mom[1]
        } else {
            tau_mis_adj <- tau_mis
        }
        if (est_syn$tau_mom[1] > 0 && tau_mom_joint[3] > 0) {
            tau_syn_adj <- tau_syn * tau_mom_joint[3] / est_syn$tau_mom[1]
        } else {
            tau_syn_adj <- tau_syn
        }
    } else {
        # Use marginal variance component tau if any group is empty
        tau_lof_adj <- tau_lof
        tau_mis_adj <- tau_mis
        tau_syn_adj <- tau_syn
    }

    # Estimate gene-level heritability
    h2_denom <- if (traitType == "binary") sum(1/v_ordered) else n_samples
    h2_lof_adj <- if (!is.null(est_lof)) compute_h2(tau_lof_adj, est_lof$tr_GtG, sigma_sq, h2_denom) else 0
    h2_mis_adj <- if (!is.null(est_mis)) compute_h2(tau_mis_adj, est_mis$tr_GtG, sigma_sq, h2_denom) else 0
    h2_syn_adj <- if (!is.null(est_syn)) compute_h2(tau_syn_adj, est_syn$tr_GtG, sigma_sq, h2_denom) else 0

    # Obtain effect size jointly if all groups are non-empty
    if ((tau_lof_adj > 0) & (tau_mis_adj > 0) & (tau_syn_adj > 0)) {
        post_beta <- calculate_joint_blup(
            G1 = lof_mat, G2 = mis_mat, G3 = syn_mat,
            tau1 = tau_lof_adj, tau2 = tau_mis_adj, tau3 = tau_syn_adj,
            psi = as.numeric(sigma_sq), y = y_vec
        )
    } else {
        # If not, calculate beta marginally
        post_beta <- as.vector(rbind(
            if (!is.null(est_lof)) est_lof$post_beta else NULL,
            if (!is.null(est_mis)) est_mis$post_beta else NULL,
            if (!is.null(est_syn)) est_syn$post_beta else NULL
        ))
    }

    post_beta <- as.vector(post_beta)

    if (isTRUE(apply_AR)) {
        tau1 = tau_lof_adj
        tau2 = tau_mis_adj
        tau3 = tau_syn_adj
        Sigma1 = rep(1, ncol(lof_mat))
        Sigma2 = rep(1, ncol(mis_mat))
        Sigma3 = rep(1, ncol(syn_mat))
        psi = as.numeric(sigma_sq)
        lambda <- psi / c(tau1 * Sigma1, tau2 * Sigma2, tau3 * Sigma3)
        result_AR <- adaptive_ridge(G_reordered, y_vec, lambda, q = 0, delta = 1e-5, gamma = 2, max_iter = 5, tol = 0.01, sigma_sq = sigma_sq)
        post_beta <- as.vector(result_AR$beta)
        se_beta <- as.vector(result_AR$se_beta)
        PEV <- se_beta^2
    } else {
        # Obtain prediction error variance (PEV)
        PEV_lof <- if (!is.null(est_lof)) compute_pev(est_lof$GtG, sigma_sq, tau_lof_adj) else NULL
        PEV_mis <- if (!is.null(est_mis)) compute_pev(est_mis$GtG, sigma_sq, tau_mis_adj) else NULL
        PEV_syn <- if (!is.null(est_syn)) compute_pev(est_syn$GtG, sigma_sq, tau_syn_adj) else NULL
        PEV <- abs(c(PEV_lof, PEV_mis, PEV_syn))
    }

    # Apply Firth bias correction for binary phenotype
    if (traitType == "binary") {
        # Align y and offset to G_reordered sample order
        firth_model_idx <- match(rownames(G_reordered), modglmm$sampleID)
        y_binary_vec <- modglmm$y[firth_model_idx]
        offset1 <- modglmm$linear.predictors[firth_model_idx] - modglmm$coefficients[1]

        out_single_wL2_lof_sparse <- numeric(0)
        out_single_wL2_mis_sparse <- numeric(0)
        out_single_wL2_syn_sparse <- numeric(0)

        # LoF
        if (lof_ncol > 0) {
            G.lof.sp <- as(G_reordered[, 1:lof_ncol, drop = F], "sparseMatrix")
            out_single_wL2_lof_sparse <- Run_Firth_MultiVar_Single(
                G.lof.sp, modglmm$obj.noK, y_binary_vec, offset1,
                lof_ncol, l2.var = 1/(2*tau_lof_adj),
                Is.Fast = FALSE, Is.Sparse = TRUE)[, 2]
        }

        # mis
        if (mis_ncol > 0) {
            mis_cols <- (lof_ncol + 1):(lof_ncol + mis_ncol)
            G.mis.sp <- as(G_reordered[, mis_cols, drop = F], "sparseMatrix")
            out_single_wL2_mis_sparse <- Run_Firth_MultiVar_Single(
                G.mis.sp, modglmm$obj.noK, y_binary_vec, offset1,
                mis_ncol, l2.var = 1/(2*tau_mis_adj),
                Is.Fast = FALSE, Is.Sparse = TRUE)[, 2]
        }

        # syn
        if (syn_ncol > 0) {
            syn_cols <- (lof_ncol + mis_ncol + 1):(lof_ncol + mis_ncol + syn_ncol)
            G.syn.sp <- as(G_reordered[, syn_cols, drop = F], "sparseMatrix")
            out_single_wL2_syn_sparse <- Run_Firth_MultiVar_Single(
                G.syn.sp, modglmm$obj.noK, y_binary_vec, offset1,
                syn_ncol, l2.var = 1/(2*tau_syn_adj),
                Is.Fast = FALSE, Is.Sparse = TRUE)[, 2]
        }

        beta_firth <- c(out_single_wL2_lof_sparse, out_single_wL2_mis_sparse, out_single_wL2_syn_sparse)
        effect <- ifelse(abs(post_beta) < log(2), post_beta, beta_firth)
    } else {
        effect <- post_beta
    }

    effect <- as.vector(effect)
    # output related to single-variant effect size
    variant <- colnames(G)

    # Build annotation and is_UR vectors
    anno_vec <- c(rep("lof", lof_ncol), rep("missense", mis_ncol), rep("synonymous", syn_ncol))
    is_UR <- grepl("_UR$", variant)

    # Parse variant IDs (chr:pos:ref:alt), UR rows get NA
    parsed <- do.call(rbind, lapply(variant, function(v) {
        parts <- strsplit(v, ":")[[1]]
        if (length(parts) >= 4) {
            data.frame(pos = as.integer(parts[2]), ref = parts[3], alt = parts[4],
                       stringsAsFactors = FALSE)
        } else {
            data.frame(pos = NA_integer_, ref = NA_character_, alt = NA_character_,
                       stringsAsFactors = FALSE)
        }
    }))

    # MAC / MAF per variant
    MAC_vec <- as.numeric(colSums(G_reordered))
    MAF_vec <- MAC_vec / (2 * n_samples)
    # UR rows: set MAC/MAF to NA (collapsed indicator, not a single variant)
    MAC_vec[is_UR] <- NA
    MAF_vec[is_UR] <- NA

    SE <- sqrt(PEV)

    # Determine sign of gene-level effect
    if (lof_ncol == 1) {
        sgn <- sign(post_beta[1])
    } else if (lof_ncol == 0) {
        message("LoF variant does not exist, so the sign of the effect size is calculated by the sign of other groups.")
        MAC_sign <- colSums(G_reordered)
        sgn <- sign(sum(MAC_sign * post_beta))
    } else {
        MAC_sign <- colSums(G_reordered[,c(1:(lof_ncol)), drop = F])
        sgn <- sign(sum(MAC_sign * post_beta[c(1:(lof_ncol))]))
    }

    # is_flipped per variant
    flipped_var_all <- rbind(lof_flipped, mis_flipped, syn_flipped)
    is_flipped_vec <- variant %in% flipped_var_all[is_flipped == TRUE, ]$var_list

    # Build effect output data.frame
    effect_out <- data.frame(
        gene = geneName,
        chr = chrom,
        pos = parsed$pos,
        ref = parsed$ref,
        alt = parsed$alt,
        variant = variant,
        annotation = anno_vec,
        is_UR = is_UR,
        MAC = MAC_vec,
        MAF = MAF_vec,
        effect = as.numeric(effect),
        SE = SE,
        PEV = PEV,
        is_flipped = is_flipped_vec,
        stringsAsFactors = FALSE
    )

    # Flip sign for allele-flipped variants
    effect_out[is_flipped_vec, "effect"] <- -effect_out[is_flipped_vec, "effect"]

    # Build h2 output data.frame (vertical format)
    h2_all <- sum(h2_lof_adj, h2_mis_adj, h2_syn_adj) * sgn
    h2_out <- data.frame(
        gene = geneName,
        chr = chrom,
        group = c("LoF", "mis", "syn", "all"),
        h2 = c(h2_lof_adj, h2_mis_adj, h2_syn_adj, abs(h2_all)),
        h2_signed = c(h2_lof_adj * sgn, h2_mis_adj * sgn, h2_syn_adj * sgn, h2_all),
        tau_adj = c(tau_lof_adj, tau_mis_adj, tau_syn_adj, NA),
        n_variants = c(lof_ncol, mis_ncol, syn_ncol, lof_ncol + mis_ncol + syn_ncol),
        stringsAsFactors = FALSE
    )

    effect_outname <- paste0(outputPrefix, "_effect.txt")
    h2_outname <- paste0(outputPrefix, "_h2.txt")

    write.table(effect_out, effect_outname, row.names = FALSE, quote = FALSE, sep = "\t")
    write.table(h2_out, h2_outname, row.names = FALSE, quote = FALSE, sep = "\t")

    message("Analysis completed.")
}
