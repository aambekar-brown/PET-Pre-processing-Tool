# ReportGenerator_UPDATED.R
suppressPackageStartupMessages({
  library(oro.nifti)
  library(tools)
})

# Simple paste0 operator (no extra deps)
"%+%" <- function(a, b) paste0(a, b)

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

.nifti_header_info <- function(img) {
  d <- dim(img)
  pix <- tryCatch(as.numeric(img@pixdim[2:4]), error = function(e) rep(NA_real_, 3))
  voxvol_mm3 <- if (all(is.finite(pix))) prod(pix) else NA_real_
  list(
    dims = paste(d, collapse = " x "),
    ndim = length(d),
    pixdim_mm = pix,
    voxvol_mm3 = voxvol_mm3,
    datatype = tryCatch(as.character(img@datatype), error = function(e) NA_character_),
    bitpix = tryCatch(as.numeric(img@bitpix), error = function(e) NA_real_),
    qform_code = tryCatch(as.numeric(img@qform_code), error = function(e) NA_real_),
    sform_code = tryCatch(as.numeric(img@sform_code), error = function(e) NA_real_)
  )
}

.nifti_intensity_summary <- function(img) {
  v_all <- as.numeric(img)
  finite <- is.finite(v_all)
  v <- v_all[finite]
  if (!length(v)) {
    return(list(
      min=NA_real_, max=NA_real_, mean=NA_real_, sd=NA_real_,
      p01=NA_real_, p05=NA_real_, p25=NA_real_, p50=NA_real_, p75=NA_real_, p95=NA_real_, p99=NA_real_,
      n_voxels=length(v_all),
      n_finite=0L, n_nan=sum(is.nan(v_all)), n_inf=sum(is.infinite(v_all)),
      frac_zero=NA_real_, frac_negative=NA_real_
    ))
  }
  qs <- as.numeric(stats::quantile(v, c(0.01,0.05,0.25,0.50,0.75,0.95,0.99), na.rm=TRUE, names=FALSE))
  list(
    min=min(v), max=max(v), mean=mean(v), sd=stats::sd(v),
    p01=qs[1], p05=qs[2], p25=qs[3], p50=qs[4], p75=qs[5], p95=qs[6], p99=qs[7],
    n_voxels=length(v_all),
    n_finite=length(v), n_nan=sum(is.nan(v_all)), n_inf=sum(is.infinite(v_all)),
    frac_zero=mean(v_all == 0, na.rm=TRUE),
    frac_negative=mean(v_all < 0 & finite, na.rm=TRUE)
  )
}

.binary_mask_from_img <- function(img, thr = 0) {
  (img > thr) & is.finite(img)
}

.mask_volume_ml <- function(mask, pixdim_mm) {
  if (length(pixdim_mm) != 3 || any(!is.finite(pixdim_mm))) return(NA_real_)
  voxvol_ml <- prod(pixdim_mm) / 1000.0
  sum(mask, na.rm=TRUE) * voxvol_ml
}

.ncc <- function(a, b, mask = NULL) {
  a <- as.numeric(a); b <- as.numeric(b)
  ok <- is.finite(a) & is.finite(b)
  if (!is.null(mask)) ok <- ok & as.logical(mask)
  a <- a[ok]; b <- b[ok]
  if (length(a) < 1000) return(NA_real_)
  if (stats::sd(a) == 0 || stats::sd(b) == 0) return(NA_real_)
  stats::cor(a, b)
}

.mutual_information <- function(a, b, mask = NULL, bins = 64) {
  a <- as.numeric(a); b <- as.numeric(b)
  ok <- is.finite(a) & is.finite(b)
  if (!is.null(mask)) ok <- ok & as.logical(mask)
  a <- a[ok]; b <- b[ok]
  if (length(a) < 1000) return(NA_real_)

  qa <- stats::quantile(a, c(0.01, 0.99), na.rm=TRUE)
  qb <- stats::quantile(b, c(0.01, 0.99), na.rm=TRUE)
  a <- pmin(pmax(a, qa[1]), qa[2])
  b <- pmin(pmax(b, qb[1]), qb[2])

  ha <- cut(a, breaks=bins, labels=FALSE, include.lowest=TRUE)
  hb <- cut(b, breaks=bins, labels=FALSE, include.lowest=TRUE)

  joint <- table(ha, hb)
  pxy <- joint / sum(joint)
  px <- rowSums(pxy)
  py <- colSums(pxy)

  nz <- pxy > 0
  sum(pxy[nz] * log(pxy[nz] / (px[row(pxy)[nz]] * py[col(pxy)[nz]])))
}

.rescale01 <- function(x, qlow = 0.01, qhigh = 0.99) {
  v <- as.numeric(x)
  v <- v[is.finite(v)]
  if (!length(v)) return(x*0)
  q <- stats::quantile(v, c(qlow, qhigh), na.rm=TRUE)
  lo <- q[1]; hi <- q[2]
  if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
    lo <- min(v); hi <- max(v)
  }
  y <- (x - lo) / (hi - lo)
  y[y < 0] <- 0
  y[y > 1] <- 1
  y[!is.finite(y)] <- 0
  y
}

.get_mid_slice <- function(img, k = NULL) {
  # Robust mid-slice extractor:
  # - works for oro.nifti objects and plain arrays
  # - supports 3D and 4D (uses first volume)
  # - never throws (returns NULL on failure)
  arr <- NULL

  if (inherits(img, "nifti")) {
    # Use the underlying array directly (more stable than subsetting nifti objects)
    arr <- tryCatch(img@.Data, error = function(e) NULL)
    if (is.null(arr)) arr <- tryCatch(as.array(img), error = function(e) NULL)
  } else {
    arr <- tryCatch(as.array(img), error = function(e) NULL)
  }

  if (is.null(arr)) return(NULL)

  d <- dim(arr)
  if (is.null(d) || length(d) < 3) return(NULL)

  # If 4D+, use the first volume
  if (length(d) >= 4 && d[4] >= 1) {
    arr <- arr[,,,1, drop = TRUE]
    d <- dim(arr)
  }

  if (is.null(d) || length(d) < 3) return(NULL)

  if (is.null(k)) k <- round(d[3] / 2)
  k <- max(1, min(d[3], k))

  sl <- arr[,,k, drop = TRUE]
  # orient for display
  t(apply(sl, 2, rev))
}



.save_mid_slice_png <- function(img, out_png, title="") {
  sl <- .get_mid_slice(img)
  if (is.null(sl)) return(NULL)
  png(out_png, width=1100, height=900, res=150)
  par(mar=c(2,2,3,2), bg="white")
  sln <- .rescale01(sl)
  image(sln, col=gray(0:255/255), axes=FALSE, asp=1, main=title)
  dev.off()
  out_png
}

.save_overlay_png <- function(bg_img, fg_img, out_png, title="", alpha = 0.50) {
  bg <- .get_mid_slice(bg_img)
  fg <- .get_mid_slice(fg_img)
  if (is.null(bg) || is.null(fg)) return(NULL)
  if (!all(dim(bg) == dim(fg))) return(NULL)

  bg01 <- .rescale01(bg)
  fg01 <- .rescale01(fg)

  r <- nrow(bg01); c <- ncol(bg01)
  pal <- grDevices::colorRampPalette(c("black","red","yellow","white"))(256)
  fg_col <- pal[pmax(1, pmin(256, floor(fg01*255)+1))]

  # col2rgb returns 3 x (r*c); reshape to 3 x r x c
  fg_rgb_mat <- grDevices::col2rgb(fg_col) / 255
  fg_rgb <- array(fg_rgb_mat, dim = c(3, r, c))

  bg_rgb <- array(rep(bg01, each=3), dim=c(3, r, c))
  out_rgb <- (1-alpha)*bg_rgb + alpha*fg_rgb
  raster <- aperm(out_rgb, c(2,3,1))

  png(out_png, width=1100, height=900, res=150)
  on.exit(dev.off(), add=TRUE)
  par(mar=c(2,2,3,2), bg="white")
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1), asp=1)
  rasterImage(raster, 0, 0, 1, 1, interpolate=FALSE)
  title(main=title)
  box()
  out_png
}



.save_mask_contour_png <- function(bg_img, mask_img, out_png, title="") {
  bg <- .get_mid_slice(bg_img)
  mk <- .get_mid_slice(mask_img)
  if (is.null(bg) || is.null(mk)) return(NULL)
  
  # Force both to 2D matrices (prevents mk being a vector with no dim)
  bg <- as.matrix(bg)
  mk <- as.matrix(mk)
  if (!all(dim(bg) == dim(mk))) return(NULL)
  
  bg01 <- .rescale01(bg)
  
  # If mask is multi-label (parcellation), draw ROI boundary edges (not mk>0)
  mk_num <- suppressWarnings(as.numeric(mk))
  mk_num[!is.finite(mk_num)] <- 0
  
  # Make sure mk_lab is a 2D matrix (same dims as bg)
  mk_lab <- matrix(
    as.integer(round(mk_num)),
    nrow = nrow(bg),
    ncol = ncol(bg),
    byrow = FALSE
  )
  
  # edge map: boundary where neighboring voxels differ
  edge <- matrix(FALSE, nrow = nrow(mk_lab), ncol = ncol(mk_lab))
  edge[-1, ] <- edge[-1, ] | (mk_lab[-1, ] != mk_lab[-nrow(mk_lab), ])
  edge[, -1] <- edge[, -1] | (mk_lab[, -1] != mk_lab[, -ncol(mk_lab)])
  edge <- edge & (mk_lab > 0)  # only draw boundaries within labeled regions
  
  png(out_png, width=1100, height=900, res=150)
  par(mar=c(2,2,3,2), bg="white")
  image(bg01, col=gray(0:255/255), axes=FALSE, asp=1, main=title)
  try(contour(edge, add=TRUE, drawlabels=FALSE, col="deepskyblue3", lwd=1), silent=TRUE)
  dev.off()
  out_png
}

.read_affine_4x4 <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  mat_txt <- readLines(path, warn=FALSE)
  vals <- suppressWarnings(as.numeric(unlist(strsplit(paste(mat_txt, collapse=" "), "\\s+"))))
  vals <- vals[is.finite(vals)]
  if (length(vals) < 16) return(NULL)
  matrix(vals[1:16], nrow=4, byrow=TRUE)
}

.affine_decompose_basic <- function(M) {
  if (is.null(M) || !is.matrix(M) || any(dim(M) != c(4,4))) return(list())
  A <- M[1:3, 1:3, drop=FALSE]
  tvec <- M[1:3, 4]
  sx <- sqrt(sum(A[,1]^2)); sy <- sqrt(sum(A[,2]^2)); sz <- sqrt(sum(A[,3]^2))
  detA <- tryCatch(det(A), error=function(e) NA_real_)
  list(
    translation_mm = tvec,
    scale_est = c(sx, sy, sz),
    detA = detA
  )
}

.file_info_row <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(list(Path=path, Exists=FALSE, SizeBytes=NA_real_, Modified=NA_character_))
  }
  fi <- file.info(path)
  list(
    Path = normalizePath(path, winslash="/", mustWork=FALSE),
    Exists = TRUE,
    SizeBytes = as.numeric(fi$size),
    Modified = as.character(fi$mtime)
  )
}

.esc <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed=TRUE)
  x <- gsub("<", "&lt;", x, fixed=TRUE)
  x <- gsub(">", "&gt;", x, fixed=TRUE)
  x <- gsub("\"", "&quot;", x, fixed=TRUE)
  x
}

.list_to_table <- function(lst, key_name = "Name", val_name = "Value") {
  if (length(lst) == 0) return("<p><em>None recorded.</em></p>")
  nms <- names(lst)
  if (is.null(nms)) nms <- rep("", length(lst))
  vals <- vapply(lst, function(v) paste0(v, collapse="; "), character(1))
  rows <- paste0("<tr><td>", .esc(nms), "</td><td>", .esc(vals), "</td></tr>", collapse="\n")
  paste0(
    "<table><thead><tr><th>", .esc(key_name), "</th><th>", .esc(val_name),
    "</th></tr></thead><tbody>", rows, "</tbody></table>"
  )
}

.fileinfo_table <- function(lst) {
  if (length(lst) == 0) return("<p><em>No files recorded.</em></p>")
  p <- unique(unlist(lst, use.names=FALSE))
  rows <- vapply(p, function(pp) {
    r <- .file_info_row(pp)
    paste0(
      "<tr><td><code>", .esc(r$Path), "</code></td>",
      "<td>", .esc(r$Exists), "</td>",
      "<td>", .esc(ifelse(is.na(r$SizeBytes), "", format(r$SizeBytes, scientific=FALSE))), "</td>",
      "<td>", .esc(ifelse(is.na(r$Modified), "", r$Modified)), "</td></tr>"
    )
  }, character(1))
  paste0(
    "<table><thead><tr><th>Path</th><th>Exists</th><th>Size (bytes)</th><th>Modified</th></tr></thead>",
    "<tbody>", paste(rows, collapse="\n"), "</tbody></table>"
  )
}

.img_to_datauri <- function(png_path) {
  if (!requireNamespace("base64enc", quietly = TRUE)) {
    install.packages("base64enc")
  }
  raw <- readBin(png_path, "raw", n = file.info(png_path)$size)
  b64 <- base64enc::base64encode(raw)
  paste0("data:image/png;base64,", b64)
}

# ------------------------------------------------------------------------------
# Main report function
# ------------------------------------------------------------------------------

generate_processing_report <- function(
    step,
    output_dir,
    input_paths = list(),
    output_paths = list(),
    parameters = list(),
    parcellation = NULL,
    logs = "",
    title = NULL,
    input_img = NULL,
    output_img = NULL,
    ref_img = NULL,
    transform_mat_path = NULL,
    extra_notes = NULL
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  if (is.null(title)) title <- paste0("TRAECR Processing Report: ", step)

  # ---------------------------
  # File info tables
  # ---------------------------
  input_files_html  <- .fileinfo_table(input_paths)
  output_files_html <- .fileinfo_table(output_paths)

  # ---------------------------
  # Image summaries
  # ---------------------------
  img_block <- function(label, img) {
    if (is.null(img) || !inherits(img, "nifti")) return(NULL)
    h <- .nifti_header_info(img)
    s <- .nifti_intensity_summary(img)
    list(
      Label = label,
      Dimensions = h$dims,
      VoxelSize_mm = paste0(sprintf("%.4g", h$pixdim_mm), collapse=" x "),
      VoxelVolume_mm3 = h$voxvol_mm3,
      DataType = h$datatype,
      Bitpix = h$bitpix,
      QForm = h$qform_code,
      SForm = h$sform_code,
      Min = s$min, Max = s$max, Mean = s$mean, SD = s$sd,
      P01 = s$p01, P05 = s$p05, P25 = s$p25, P50 = s$p50, P75 = s$p75, P95 = s$p95, P99 = s$p99,
      FiniteVoxels = s$n_finite,
      NaN_Count = s$n_nan,
      Inf_Count = s$n_inf,
      FractionZero = round(s$frac_zero, 6),
      FractionNegative = round(s$frac_negative, 6)
    )
  }

  img_summaries <- list()
  ib <- img_block("Input", input_img);   if (!is.null(ib)) img_summaries <- c(img_summaries, list(ib))
  ob <- img_block("Output", output_img); if (!is.null(ob)) img_summaries <- c(img_summaries, list(ob))
  rb <- img_block("Reference", ref_img); if (!is.null(rb)) img_summaries <- c(img_summaries, list(rb))

  img_summary_html <- if (!length(img_summaries)) {
    "<p><em>No NIfTI images were provided to the report function.</em></p>"
  } else {
    paste0(
      vapply(img_summaries, function(b) {
        "<h3>" %+% .esc(b$Label) %+% "</h3>" %+% .list_to_table(b, "Field", "Value")
      }, character(1)),
      collapse = "\n"
    )
  }

  # ---------------------------
  # QC metrics (actionable)
  # ---------------------------
  qc_metrics <- list(
    Step = step,
    Timestamp = ts,
    OutputFolder = normalizePath(output_dir, winslash="/", mustWork=FALSE)
  )

  # Empty output detection + volume
  if (!is.null(output_img) && inherits(output_img, "nifti")) {
    h <- .nifti_header_info(output_img)
    out_mask <- .binary_mask_from_img(output_img, thr=0)
    qc_metrics$OutputNonZeroVoxels  <- sum(out_mask, na.rm=TRUE)
    qc_metrics$OutputNonZeroPercent <- round(100 * mean(out_mask, na.rm=TRUE), 3)
    qc_metrics$OutputNonZeroVolume_mL <- .mask_volume_ml(out_mask, h$pixdim_mm)
    qc_metrics$SuccessFlag <- if (sum(out_mask, na.rm=TRUE) > 5000) "PASS" else "FAIL/EMPTY"
  }

  # Similarity: input vs output (only if same grid)
  if (!is.null(input_img) && inherits(input_img, "nifti") &&
      !is.null(output_img) && inherits(output_img, "nifti") &&
      all(dim(input_img)[1:3] == dim(output_img)[1:3])) {
    m <- .binary_mask_from_img(output_img, 0)
    qc_metrics$NCC_Input_vs_Output_in_OutputMask <- .ncc(input_img, output_img, m)
    qc_metrics$MI_Input_vs_Output_in_OutputMask  <- .mutual_information(input_img, output_img, m, bins=64)
  }

  # Similarity: output vs reference (only if same grid)
  if (!is.null(output_img) && inherits(output_img, "nifti") &&
      !is.null(ref_img) && inherits(ref_img, "nifti") &&
      all(dim(output_img)[1:3] == dim(ref_img)[1:3])) {
    mask <- .binary_mask_from_img(output_img, 0) & .binary_mask_from_img(ref_img, 0)
    qc_metrics$NCC_Output_vs_Reference <- .ncc(output_img, ref_img, mask)
    qc_metrics$MI_Output_vs_Reference  <- .mutual_information(output_img, ref_img, mask, bins=64)
  }

  # Affine (registration/coreg)
  if (!is.null(transform_mat_path) && file.exists(transform_mat_path)) {
    M <- .read_affine_4x4(transform_mat_path)
    qc_metrics$TransformMatrixPath <- normalizePath(transform_mat_path, winslash="/", mustWork=FALSE)
    if (!is.null(M)) {
      de <- .affine_decompose_basic(M)
      qc_metrics$AffineDeterminant_3x3 <- de$detA
      if (!is.null(de$translation_mm)) qc_metrics$Translation_mm <- paste0(sprintf("%.4g", de$translation_mm), collapse=", ")
      if (!is.null(de$scale_est)) qc_metrics$ScaleEstimates <- paste0(sprintf("%.4g", de$scale_est), collapse=", ")
    }
  }

  # Brain extraction: estimate brain volume from output (non-zero)
  if (grepl("brain extraction", step, ignore.case=TRUE) &&
      !is.null(output_img) && inherits(output_img, "nifti")) {
    h <- .nifti_header_info(output_img)
    m <- .binary_mask_from_img(output_img, 0)
    qc_metrics$EstimatedBrainVolume_mL <- .mask_volume_ml(m, h$pixdim_mm)
  }

  # Harmonization/normalization: distribution shift summary
  if (grepl("combat|ravel|harmon", step, ignore.case=TRUE) &&
      !is.null(input_img) && inherits(input_img, "nifti") &&
      !is.null(output_img) && inherits(output_img, "nifti") &&
      all(dim(input_img)[1:3] == dim(output_img)[1:3])) {
    m <- .binary_mask_from_img(output_img, 0)
    a <- as.numeric(input_img); b <- as.numeric(output_img)
    ok <- is.finite(a) & is.finite(b) & as.logical(m)
    a <- a[ok]; b <- b[ok]
    if (length(a) > 1000) {
      qc_metrics$DeltaMean_inMask <- mean(b) - mean(a)
      qc_metrics$DeltaSD_inMask   <- stats::sd(b) - stats::sd(a)
      qa <- stats::quantile(a, c(0.05,0.50,0.95), na.rm=TRUE)
      qb <- stats::quantile(b, c(0.05,0.50,0.95), na.rm=TRUE)
      qc_metrics$DeltaP05_inMask <- qb[1] - qa[1]
      qc_metrics$DeltaP50_inMask <- qb[2] - qa[2]
      qc_metrics$DeltaP95_inMask <- qb[3] - qa[3]
    }
  }

  # ---------------------------
  # Optional parcellation (ONLY for template registration + MRI–PET co-registration)
  # ---------------------------
  parc_section_html <- ""
  step_allows_parc <- (grepl("template registration", step, ignore.case = TRUE) ||
                       (grepl("mri", step, ignore.case = TRUE) && grepl("pet", step, ignore.case = TRUE) &&
                        grepl("co[-–]?registration", step, ignore.case = TRUE)))

  # Prepare a small metadata block (no effect unless parcellation was passed in)
  atlas_path <- NA_character_
  csv_path   <- NA_character_
  note_txt   <- ""
  parc_performed <- FALSE
  parc_ref_path <- NA_character_

  if (step_allows_parc && !is.null(parcellation) && is.list(parcellation)) {
    
    parc_ref_path <- tryCatch(as.character(parcellation$ref_path), error = function(e) NA_character_)
    
    # Prefer PET/reference-space products if available
    atlas_labels_path <- tryCatch(as.character(parcellation$atlas_labels_path), error = function(e) NA_character_)
    resampled_atlas_path <- tryCatch(as.character(parcellation$resampled_atlas_path), error = function(e) NA_character_)
    atlas_orig_path <- tryCatch(as.character(parcellation$atlas_path), error = function(e) NA_character_)
    
    # Choose best existing file for visualization (labelmap in PET/ref space is best)
    viz_atlas_path <- NA_character_
    if (is.character(atlas_labels_path) && nzchar(atlas_labels_path) && file.exists(atlas_labels_path)) {
      viz_atlas_path <- atlas_labels_path
    } else if (is.character(resampled_atlas_path) && nzchar(resampled_atlas_path) && file.exists(resampled_atlas_path)) {
      viz_atlas_path <- resampled_atlas_path
    } else if (is.character(atlas_orig_path) && nzchar(atlas_orig_path) && file.exists(atlas_orig_path)) {
      viz_atlas_path <- atlas_orig_path
    }
    
    csv_path   <- tryCatch(as.character(parcellation$csv_path), error = function(e) NA_character_)
    note_txt   <- tryCatch(as.character(parcellation$note),     error = function(e) "")
    parc_performed <- tryCatch(isTRUE(parcellation$performed),  error = function(e) FALSE)
    roi_dir <- tryCatch(as.character(parcellation$roi_dir),     error = function(e) NA_character_)
    
    # Store chosen atlas path in the existing variable name used downstream
    atlas_path <- viz_atlas_path
    
    parc_meta <- list(
      performed = parc_performed,
      atlas_labels_path = if (!is.na(atlas_labels_path)) atlas_labels_path else "",
      resampled_atlas_path = if (!is.na(resampled_atlas_path)) resampled_atlas_path else "",
      atlas_used_for_viz = if (!is.na(viz_atlas_path)) viz_atlas_path else "",
      roi_dir = if (!is.na(roi_dir)) roi_dir else "",
      csv_path = if (!is.na(csv_path)) csv_path else "",
      note = if (!is.na(note_txt)) note_txt else ""
    )
    
    parc_section_html <- "<h2>Parcellation</h2>" %+% .list_to_table(parc_meta, "Name", "Value")
  }


  # ---------------------------
  # Figures
  # ---------------------------
  figures <- character(0)

  # Parcellation visualization (ONLY for template registration + MRI–PET co-registration)
  if (step_allows_parc && isTRUE(parc_performed) &&
      is.character(atlas_path) && nzchar(atlas_path) && file.exists(atlas_path)) {
    
    atlas_img <- tryCatch(readNIfTI(atlas_path, reorient = FALSE), error = function(e) NULL)
    
    # Prefer to overlay parcellation on its own reference image (PET space)
    bg_img <- NULL
    if (is.character(parc_ref_path) && nzchar(parc_ref_path) && file.exists(parc_ref_path)) {
      bg_img <- tryCatch(readNIfTI(parc_ref_path, reorient = FALSE), error = function(e) NULL)
    }
    
    # Fallback to output_img if no ref image available
    if (is.null(bg_img) && !is.null(output_img) && inherits(output_img, "nifti")) {
      bg_img <- output_img
    }
    
    if (!is.null(atlas_img) && inherits(atlas_img, "nifti") &&
        !is.null(bg_img) && inherits(bg_img, "nifti") &&
        all(dim(atlas_img)[1:3] == dim(bg_img)[1:3])) {
      
      figures <- c(figures, .save_mask_contour_png(
        bg_img, atlas_img,
        file.path(output_dir, "report_parcellation_contour_on_petspace.png"),
        "Parcellation contour on PET (mid-slice)"
      ))
    }
  }

  if (!is.null(input_img) && inherits(input_img, "nifti")) {
    figures <- c(figures, .save_mid_slice_png(input_img, file.path(output_dir, "report_input_mid_slice.png"), "Input (mid-slice)"))
  }
  if (!is.null(output_img) && inherits(output_img, "nifti")) {
    figures <- c(figures, .save_mid_slice_png(output_img, file.path(output_dir, "report_output_mid_slice.png"), "Output (mid-slice)"))
  }
  if (!is.null(ref_img) && inherits(ref_img, "nifti")) {
    figures <- c(figures, .save_mid_slice_png(ref_img, file.path(output_dir, "report_reference_mid_slice.png"), "Reference (mid-slice)"))
  }
  if (!is.null(input_img) && inherits(input_img, "nifti") &&
      !is.null(output_img) && inherits(output_img, "nifti") &&
      all(dim(input_img)[1:3] == dim(output_img)[1:3])) {
    figures <- c(figures, .save_overlay_png(input_img, output_img, file.path(output_dir, "report_overlay_input_bg_output_fg.png"),
                                           "Overlay: input (gray) + output (heat)", alpha=0.45))
  }
  if (grepl("brain extraction", step, ignore.case=TRUE) &&
      !is.null(input_img) && inherits(input_img, "nifti") &&
      !is.null(output_img) && inherits(output_img, "nifti") &&
      all(dim(input_img)[1:3] == dim(output_img)[1:3])) {
    figures <- c(figures, .save_mask_contour_png(input_img, output_img, file.path(output_dir, "report_brainmask_contour_on_input.png"),
                                                 "Brain mask contour on input (mid-slice)"))
  }

  figures <- Filter(function(x) !is.null(x) && file.exists(x), figures)
  figs_html <- if (!length(figures)) {
    "<p><em>No figures generated.</em></p>"
  } else {
    paste0(
      vapply(figures, function(f) {
        uri <- .img_to_datauri(f)
        "<div class='fig'><div class='figcap'>" %+% .esc(basename(f)) %+% "</div><img src='" %+% uri %+% "' /></div>"
      }, character(1)),
      collapse = "\n"
    )
  }

  # ---------------------------
  # Render self-contained HTML
  # ---------------------------
  css <- "
    body { font-family: -apple-system,BlinkMacSystemFont,Segoe UI,Roboto,Arial,sans-serif; margin: 24px; }
    h1 { margin: 0 0 6px 0; }
    .meta { color:#555; margin-bottom: 18px; }
    h2 { margin-top: 22px; }
    h3 { margin-top: 14px; }
    table { border-collapse: collapse; width: 100%; margin: 10px 0 18px 0; }
    th, td { border: 1px solid #ddd; padding: 8px; vertical-align: top; }
    th { background: #f6f6f6; text-align: left; }
    code { background:#f3f3f3; padding:2px 4px; border-radius:4px; }
    pre { background:#0b0b0b; color:#eaeaea; padding:12px; border-radius:8px; overflow:auto; }
    .fig img { max-width: 100%; border: 1px solid #ddd; border-radius: 6px; }
    .figcap { color:#444; font-size: 13px; margin: 8px 0; }
    .note { padding: 10px 12px; background: #fff7d6; border: 1px solid #f0d58c; border-radius: 8px; }
  "

  html <- "<!doctype html><html><head><meta charset='utf-8'>" %+%
    "<title>" %+% .esc(title) %+% "</title>" %+%
    "<style>" %+% css %+% "</style></head><body>" %+%
    "<h1>" %+% .esc(title) %+% "</h1>" %+%
    "<div class='meta'><b>Step:</b> " %+% .esc(step) %+% " &nbsp; | &nbsp; <b>Timestamp:</b> " %+% .esc(ts) %+% "</div>" %+%
    "<h2>What this report covers</h2>" %+%
    "<div class='note'>Actionable QC for an image pre-processing step:<ul><li>Header/geometry checks</li><li>Intensity sanity checks (NaNs/Infs/negatives)</li><li>Empty-output detection</li><li>Similarity metrics (NCC/MI when grids match)</li><li>Registration affine sanity (translation/scale/determinant)</li><li>Quick mid-slice + overlay figures</li></ul></div>" %+%
    "<h2>Inputs</h2>" %+% input_files_html %+%
    "<h2>Outputs</h2>" %+% output_files_html %+%
    "<h2>Parameters</h2>" %+% .list_to_table(parameters, "Name", "Value") %+%
    "<h2>Image summaries</h2>" %+% img_summary_html %+%
    "<h2>QC metrics</h2>" %+% .list_to_table(qc_metrics, "Metric", "Value") %+%
    parc_section_html %+%
    "<h2>Quick QC figures</h2>" %+% figs_html %+%
    (if (!is.null(extra_notes) && nzchar(extra_notes)) "<h2>Notes</h2><pre>" %+% .esc(extra_notes) %+% "</pre>" else "") %+%
    "<h2>Log</h2><pre>" %+% .esc(logs) %+% "</pre>" %+%
    "<h2>Session Info</h2><pre>" %+% .esc(paste(capture.output(sessionInfo()), collapse="\n")) %+% "</pre>" %+%
    "</body></html>"

  out_html <- file.path(output_dir, "processing_report.html")
  writeLines(html, out_html, useBytes = TRUE)
  return(out_html)
}

# ------------------------------------------------------------------------------
# Batch runner (optional helper)
# ------------------------------------------------------------------------------

generate_reports_for_existing_outputs <- function(base_dir, step_name = "unknown") {
  dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
  out <- character(0)
  for (d in dirs) {
    out_html <- tryCatch({
      generate_processing_report(
        step = step_name,
        output_dir = d,
        input_paths = list(),
        output_paths = list(),
        parameters = list(),
        logs = "Batch regeneration (no images/metrics loaded)."
      )
    }, error = function(e) NA_character_)
    out <- c(out, out_html)
  }
  out
}
