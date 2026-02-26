###############################################################################
# KIC CALCIUM ANALYSIS PIPELINE — STEP 3
# -----------------------------------------------------------------------------
#
# PURPOSE:
#   Interactive GUI-based plotting and statistical analysis tool
#   for calcium datasets generated in Step 2a (single-batch)
#   or Step 2b (multi-batch merge).
#
#   The script:
#     1. Prompts the user to select an input merged CSV file
#        and an output folder (tcltk GUI)
#     2. Allows optional selection of specific batches and/or groups
#     3. Allows user-defined group ordering (↑ / ↓ reordering in GUI)
#     4. Applies strict Peaks quality control:
#          - Uses fixed Num.Peaks column
#          - Validates expected peak count based on manual Hz × 10 sec window
#          - Applies user-defined tolerance threshold
#     5. Performs automatic numeric coercion and column harmonization
#     6. Applies within-group outlier detection (1.5×IQR rule)
#     7. Generates publication-ready bar plots:
#          - Mean ± SD bars
#          - Individual well dots (color-coded by Batch)
#     8. Performs statistical testing:
#          - Welch t-test (2 groups)
#          - One-way ANOVA + Tukey post-hoc (≥3 groups)
#          - Significance annotation (*, **, ***, ****)
#     9. Exports each metric as both PNG and PDF (vector format)
#    10. Generates diagnostic and quality-control reports
#
# OUTPUT:
#   An output folder containing:
#     *.png                      (bar plots with dots and stats)
#     *.pdf                      (vector versions of all plots)
#     plot_diagnostics.csv       (per-metric plotting summary)
#     outliers_removed.csv       (IQR-filtered data points)
#     peak_filter_log.csv        (Peaks QC exclusion log)
#
# NOTES:
#   - Uses tcltk for interactive GUI selection (macOS compatible).
#   - Enforces user-defined group order in plots.
#   - Designed for interactive exploration and final figure generation.
#   - Intended as Step 3 of the KIC calcium processing pipeline.
#   - Expects as input a merged dataset produced by:
#         • Step 2a (single-batch merged.csv), or
#         • Step 2b (merged_multibatch.csv).
#
# AUTHORS:
#   Michele Buono
#   Talitha Spanjersberg
#   Nikki Scheen
#   Nina van der Wilt
#   Regenerative Medicine Center Utrecht (2026)
###############################################################################

suppressPackageStartupMessages({
  library(tcltk)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(tibble)
})

# --------------------------- Fixed defaults -----------------------------------

PEAKS_COL    <- "Num.Peaks"
WINDOW_SEC   <- 10
OUTLIER_COEF <- 1.5

# --------------------------- Helpers -----------------------------------------

safe_message <- function(...) cat(paste0(..., "\n"))

choose_file_gui <- function(caption = "Select CSV file") {
  f <- tclvalue(tkgetOpenFile(title = caption, filetypes = "{{CSV Files} {.csv}} {{All files} {*}}"))
  if (is.na(f) || !nzchar(f)) return(NA_character_)
  normalizePath(f, winslash = "/", mustWork = TRUE)
}

choose_dir_gui <- function(caption = "Select output folder") {
  d <- tclvalue(tkchooseDirectory(title = caption))
  if (is.na(d) || !nzchar(d)) return(NA_character_)
  normalizePath(d, winslash = "/", mustWork = TRUE)
}

sanitize_filename <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)

is_integer64 <- function(x) inherits(x, "integer64")

coerce_numeric_like <- function(x){
  if (is.numeric(x))  return(x)
  if (is_integer64(x)) return(as.numeric(x))
  if (is.factor(x))   x <- as.character(x)
  if (!is.character(x)) return(x)
  
  x <- trimws(gsub("\u00A0","",x, fixed = TRUE))
  x[x == ""] <- NA
  
  base <- gsub("%", "", x, fixed = TRUE)
  A <- gsub(",", ".", base, fixed = TRUE); A <- gsub("[^0-9eE+\\-.]", "", A, perl = TRUE)
  B <- gsub(",", "", base, fixed = TRUE);  B <- gsub("[^0-9eE+\\-.]", "", B, perl = TRUE)
  
  xA <- suppressWarnings(as.numeric(A))
  xB <- suppressWarnings(as.numeric(B))
  
  if (sum(!is.na(xA)) >= sum(!is.na(xB))) xA else xB
}

flag_iqr <- function(x, coef = 1.5){
  q <- stats::quantile(x, probs = c(0.25, 0.75), na.rm = TRUE, type = 7)
  iqr <- q[[2]] - q[[1]]
  if (!is.finite(iqr) || iqr == 0) return(rep(FALSE, length(x)))
  lo <- q[[1]] - coef*iqr
  hi <- q[[2]] + coef*iqr
  x < lo | x > hi
}

sig_label <- function(p){
  if (is.na(p)) return("ns")
  if (p < 1e-4) "****" else if (p < 1e-3) "***" else if (p < 1e-2) "**" else if (p < 0.05) "*" else "ns"
}

palette_choices <- c("Pastel", "Vibrant", "Dark", "Okabe-Ito")

get_palette <- function(n, choice){
  if (n <= 0) return(character())
  choice <- as.character(choice)
  if (choice == "Okabe-Ito") {
    ok <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#000000")
    return(rep(ok, length.out = n))
  }
  pal <- switch(choice,
                "Pastel"  = "Pastel 1",
                "Vibrant" = "Set 2",
                "Dark"    = "Dark 2",
                "Set 3")
  if (exists("hcl.colors", mode="function")) {
    return(grDevices::hcl.colors(n, palette = pal))
  }
  grDevices::rainbow(n)
}

read_colnames_only <- function(path){
  x <- suppressMessages(readr::read_csv(path, n_max = 0, show_col_types = FALSE))
  names(x)
}

read_unique_levels_fast <- function(path, col){
  tryCatch({
    bx <- suppressMessages(readr::read_csv(
      path,
      col_select = all_of(col),
      show_col_types = FALSE,
      progress = FALSE
    ))
    sort(unique(as.character(bx[[col]])))
  }, error = function(e) NULL)
}

get_listbox_selection <- function(lb, values_vec){
  sel <- as.integer(tkcurselection(lb))
  if (!length(sel)) return(character())
  values_vec[sel + 1]
}

# --- NEW: listbox order + move up/down (for Group order control) --------------
lb_get_items <- function(lb){
  n <- as.integer(tkindex(lb, "end"))
  if (!is.finite(n) || n <= 0) return(character())
  sapply(0:(n-1), function(i) as.character(tkget(lb, i)))
}

lb_move_selected <- function(lb, direction = c("up","down")){
  direction <- match.arg(direction)
  sel <- as.integer(tkcurselection(lb))
  if (!length(sel)) return(invisible(NULL))
  
  items <- lb_get_items(lb)
  if (!length(items)) return(invisible(NULL))
  
  if (direction == "up") {
    for (i in sel) {
      if (i <= 0) next
      tmp <- items[i+1]; items[i+1] <- items[i]; items[i] <- tmp
    }
    new_sel <- pmax(sel - 1, 0)
  } else {
    for (i in rev(sel)) {
      if (i >= length(items)-1) next
      tmp <- items[i+1]; items[i+1] <- items[i+2]; items[i+2] <- tmp
    }
    new_sel <- pmin(sel + 1, length(items)-1)
  }
  
  tkdelete(lb, 0, "end")
  for (x in items) tkinsert(lb, "end", x)
  
  tkselection.clear(lb, 0, "end")
  for (i in new_sel) tkselection.set(lb, i)
  
  invisible(NULL)
}

get_listbox_selection_in_current_order <- function(lb){
  items <- lb_get_items(lb)
  sel <- as.integer(tkcurselection(lb))
  if (!length(sel)) return(character())
  items[sel + 1]
}

save_plot_png_pdf <- function(p, base_path_no_ext, width=7, height=4.6, dpi=300){
  # PNG
  ggsave(
    filename = paste0(base_path_no_ext, ".png"),
    plot = p,
    width = width,
    height = height,
    dpi = dpi
  )
  
  # PDF (vector). Try cairo_pdf; fallback to regular pdf device if needed.
  ok <- TRUE
  tryCatch({
    ggsave(
      filename = paste0(base_path_no_ext, ".pdf"),
      plot = p,
      width = width,
      height = height,
      device = cairo_pdf
    )
  }, error = function(e){
    ok <<- FALSE
  })
  if (!ok) {
    ggsave(
      filename = paste0(base_path_no_ext, ".pdf"),
      plot = p,
      width = width,
      height = height,
      device = "pdf"
    )
  }
}

# --------------------------- GUI ---------------------------------------------

tt <- tktoplevel()
tkwm.title(tt, "Script 3 — Bar+SD plots (GUI)")

v_infile   <- tclVar("")
v_outdir   <- tclVar("")
v_outname  <- tclVar("plots_barstats_SD_all")

v_groupcol <- tclVar("Group")
v_batchcol <- tclVar("Batch")

v_batch_mode <- tclVar("all")   # all / select
v_group_mode <- tclVar("all")   # all / select

v_hz      <- tclVar("1")
v_tol_pct <- tclVar("30")

v_palette <- tclVar("Pastel")
v_status  <- tclVar("Ready.")

row <- 0

tkgrid(tklabel(tt, text="1) Input CSV:"), row=row, column=0, sticky="w", padx=6, pady=4)
tkgrid(tkentry(tt, textvariable=v_infile, width=55), row=row, column=1, sticky="we", padx=6, pady=4)

batch_values <- character()
group_values <- character()
lb_batch <- NULL
lb_group <- NULL

populate_listboxes <- function(file){
  nms <- read_colnames_only(file)
  bcol <- tclvalue(v_batchcol)
  gcol <- tclvalue(v_groupcol)
  
  if (bcol %in% nms) {
    batch_values <<- read_unique_levels_fast(file, bcol)
    if (is.null(batch_values) || !length(batch_values)) batch_values <<- character()
  } else {
    batch_values <<- "Batch1"
  }
  
  if (gcol %in% nms) {
    group_values <<- read_unique_levels_fast(file, gcol)
    if (is.null(group_values) || !length(group_values)) group_values <<- character()
  } else {
    group_values <<- character()
  }
  
  tkdelete(lb_batch, 0, "end")
  for (x in batch_values) tkinsert(lb_batch, "end", x)
  if (length(batch_values)) tkselection.set(lb_batch, 0)
  
  # IMPORTANT: listbox order becomes the "reorderable" donor order.
  tkdelete(lb_group, 0, "end")
  for (x in group_values) tkinsert(lb_group, "end", x)
  if (length(group_values)) tkselection.set(lb_group, 0)
}

tkgrid(tkbutton(tt, text="Browse...", command=function(){
  f <- choose_file_gui("Select the filtered CSV from Script 2a or 2b")
  if (is.na(f)) return(invisible(NULL))
  tclvalue(v_infile) <- f
  populate_listboxes(f)
}), row=row, column=2, padx=6, pady=4)
row <- row + 1

tkgrid(tklabel(tt, text="2) Output folder (parent):"), row=row, column=0, sticky="w", padx=6, pady=4)
tkgrid(tkentry(tt, textvariable=v_outdir, width=55), row=row, column=1, sticky="we", padx=6, pady=4)
tkgrid(tkbutton(tt, text="Choose...", command=function(){
  d <- choose_dir_gui("Choose where to create the plot output folder")
  if (!is.na(d)) tclvalue(v_outdir) <- d
}), row=row, column=2, padx=6, pady=4)
row <- row + 1

tkgrid(tklabel(tt, text="3) Output folder name:"), row=row, column=0, sticky="w", padx=6, pady=4)
tkgrid(tkentry(tt, textvariable=v_outname, width=55), row=row, column=1, sticky="we", padx=6, pady=4)
row <- row + 1

tkgrid(tklabel(tt, text="4) Columns (normally leave defaults):"), row=row, column=0, sticky="w", padx=6, pady=4)
frm_cols <- tkframe(tt)
tkgrid(frm_cols, row=row, column=1, sticky="w")
tkgrid(tklabel(frm_cols, text="Group col:"), row=0, column=0, padx=4)
tkgrid(tkentry(frm_cols, textvariable=v_groupcol, width=12), row=0, column=1, padx=4)
tkgrid(tklabel(frm_cols, text="Batch col:"), row=0, column=2, padx=4)
tkgrid(tkentry(frm_cols, textvariable=v_batchcol, width=12), row=0, column=3, padx=4)
row <- row + 1

tkgrid(tklabel(tt, text="5) Select batches to include:"), row=row, column=0, sticky="nw", padx=6, pady=4)
frm_batch <- tkframe(tt)
tkgrid(frm_batch, row=row, column=1, sticky="w", padx=6, pady=4)
tkgrid(tkradiobutton(frm_batch, text="All batches", value="all", variable=v_batch_mode),
       row=0, column=0, sticky="w", padx=4)
tkgrid(tkradiobutton(frm_batch, text="Select batches (multi-select)", value="select", variable=v_batch_mode),
       row=1, column=0, sticky="w", padx=4)
lb_batch <- tklistbox(frm_batch, height=6, selectmode="multiple", exportselection=FALSE, width=32)
tkgrid(lb_batch, row=2, column=0, sticky="w", padx=4, pady=2)
tkgrid(tklabel(frm_batch, text="(List fills after choosing file)"), row=3, column=0, sticky="w", padx=4)
row <- row + 1

tkgrid(tklabel(tt, text="6) Select groups to include (reorder with ↑/↓):"), row=row, column=0, sticky="nw", padx=6, pady=4)
frm_group <- tkframe(tt)
tkgrid(frm_group, row=row, column=1, sticky="w", padx=6, pady=4)
tkgrid(tkradiobutton(frm_group, text="All groups", value="all", variable=v_group_mode),
       row=0, column=0, sticky="w", padx=4)
tkgrid(tkradiobutton(frm_group, text="Select groups (multi-select)", value="select", variable=v_group_mode),
       row=1, column=0, sticky="w", padx=4)

# listbox + reorder buttons side-by-side
frm_group_list <- tkframe(frm_group)
tkgrid(frm_group_list, row=2, column=0, sticky="w", padx=4, pady=2)

lb_group <- tklistbox(frm_group_list, height=6, selectmode="multiple", exportselection=FALSE, width=32)
tkgrid(lb_group, row=0, column=0, sticky="w")

frm_group_btns <- tkframe(frm_group_list)
tkgrid(frm_group_btns, row=0, column=1, sticky="n", padx=6)

tkgrid(tkbutton(frm_group_btns, text="↑", width=3,
                command=function(){ lb_move_selected(lb_group, "up") }),
       row=0, column=0, pady=2)
tkgrid(tkbutton(frm_group_btns, text="↓", width=3,
                command=function(){ lb_move_selected(lb_group, "down") }),
       row=1, column=0, pady=2)

tkgrid(tklabel(frm_group, text="(List fills after choosing file)"), row=3, column=0, sticky="w", padx=4)
row <- row + 1

tkgrid(tklabel(tt, text="7) Peaks QC (always applied):"), row=row, column=0, sticky="w", padx=6, pady=4)
frm_pk <- tkframe(tt)
tkgrid(frm_pk, row=row, column=1, sticky="w")
tkgrid(tklabel(frm_pk, text=paste0("Uses ", PEAKS_COL, "; window=", WINDOW_SEC, "s")), row=0, column=0, padx=4, sticky="w")
tkgrid(tklabel(frm_pk, text="Hz:"), row=0, column=1, padx=4, sticky="w")
tkgrid(tkentry(frm_pk, textvariable=v_hz, width=6), row=0, column=2, padx=4)
tkgrid(tklabel(frm_pk, text="Tolerance (%):"), row=0, column=3, padx=4, sticky="w")
tkgrid(tkentry(frm_pk, textvariable=v_tol_pct, width=6), row=0, column=4, padx=4)
tkgrid(tklabel(frm_pk, text="Example: 1Hz → 10 peaks; 30% → 7..13"), row=1, column=0, columnspan=6, padx=4, sticky="w")
row <- row + 1

tkgrid(tklabel(tt, text="8) Palette (group bar fills):"), row=row, column=0, sticky="w", padx=6, pady=4)
pal_widget <- if (exists("ttkcombobox")) {
  ttkcombobox(tt, textvariable = v_palette, values = palette_choices, width = 18, state = "readonly")
} else {
  tkentry(tt, textvariable = v_palette, width = 20)
}
tkgrid(pal_widget, row=row, column=1, sticky="w", padx=6, pady=4)
row <- row + 1

tkgrid(tklabel(tt, textvariable=v_status, foreground="blue"), row=row, column=1, sticky="w", padx=6, pady=4)
tkgrid.columnconfigure(tt, 1, weight=1)

# ------------------------ Plotting -------------------------------------------

make_plot_one <- function(df, var, palette_choice){
  d0 <- df |>
    transmute(Well=.data[["Well"]], Group=.data[["Group"]], Batch=.data[["Batch"]], value=.data[[var]]) |>
    tidyr::drop_na(value)
  if (!nrow(d0)) return(list(error="all NA (no finite values)"))
  
  d0 <- d0 %>% group_by(Group) %>% mutate(.rm = flag_iqr(value, OUTLIER_COEF)) %>% ungroup()
  removed <- d0 %>% filter(.rm) %>% select(Well, Group, Batch, value) %>% mutate(variable=var)
  d <- d0 %>% filter(!.rm) %>% select(-.rm)
  if (!nrow(d)) return(list(error="all removed as outliers", removed=removed))
  
  # ensure factors (use df levels)
  d$Group <- factor(d$Group, levels=levels(df$Group))
  d$Batch <- factor(d$Batch, levels=levels(df$Batch))
  
  sum_df <- d %>%
    group_by(Group) %>%
    summarise(n=dplyr::n(), mean=mean(value), sd=ifelse(n()>1, sd(value), 0), .groups="drop")
  if (!nrow(sum_df) || all(!is.finite(sum_df$mean))) return(list(error="no finite group means", removed=removed))
  
  group_levels <- levels(sum_df$Group)
  fills <- get_palette(length(group_levels), palette_choice); names(fills) <- group_levels
  batch_levels <- levels(d$Batch)
  batch_cols <- get_palette(length(batch_levels), "Vibrant"); names(batch_cols) <- batch_levels
  
  p <- ggplot(sum_df, aes(x=Group, y=mean, fill=Group)) +
    geom_col(width=0.6) +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.22) +
    geom_jitter(data=d, aes(x=Group, y=value, color=Batch),
                width=0.10, height=0, alpha=0.75, size=1.9, inherit.aes=FALSE) +
    scale_fill_manual(values=fills) +
    scale_color_manual(values=batch_cols) +
    labs(title=var, x=NULL, y=paste0(var, " (mean ± SD)"), fill="Group", color="Batch") +
    theme_classic(base_size=12)
  
  ymax_bar <- max(sum_df$mean + sum_df$sd, na.rm=TRUE)
  ymin_bar <- min(sum_df$mean - sum_df$sd, na.rm=TRUE)
  n_groups <- dplyr::n_distinct(d$Group)
  headroom <- if (n_groups==2) 0.20 else 0.35
  ymax_plot <- ymax_bar * (1 + headroom)
  ymin_plot <- min(0, ymin_bar)
  p <- p + coord_cartesian(ylim=c(ymin_plot, ymax_plot))
  
  note <- paste0("outliers removed: ", nrow(removed), " (IQR ", OUTLIER_COEF, "×)")
  
  # Stats + brackets
  ann <- NULL
  xlev <- levels(sum_df$Group)
  xmap <- setNames(seq_along(xlev), xlev)
  
  if (n_groups == 2) {
    gcounts <- table(d$Group)
    if (all(gcounts >= 2)) {
      tt <- t.test(value ~ Group, data=d, var.equal=FALSE)
      lab <- sig_label(tt$p.value)
      if (lab != "ns") ann <- rbind(ann, data.frame(x1=1,x2=2,y=ymax_bar*1.10,label=lab))
      note <- paste0(note, " | Welch t-test p=", signif(tt$p.value, 3))
    } else {
      note <- paste0(note, " | stats skipped (need ≥2 per group)")
    }
  } else if (n_groups >= 3) {
    gcounts <- table(d$Group)
    if (sum(gcounts >= 2) >= 2) {
      aovm <- aov(value ~ Group, data=d)
      anov_p <- summary(aovm)[[1]][["Pr(>F)"]][1]
      note <- paste0(note, " | ANOVA p=", signif(anov_p, 3))
      
      tk <- TukeyHSD(aovm)$Group
      tk_df <- as.data.frame(tk); tk_df$pair <- rownames(tk_df)
      sig_pairs <- tk_df %>%
        mutate(label = vapply(`p adj`, sig_label, character(1))) %>%
        filter(label != "ns")
      
      if (nrow(sig_pairs)) {
        parse_pair <- function(s) strsplit(s, "-")[[1]]
        height <- ymax_bar * 1.06
        step <- 0.06 * ymax_bar
        i <- 0L
        for (k in seq_len(nrow(sig_pairs))) {
          ab <- parse_pair(sig_pairs$pair[k])
          ann <- rbind(ann, data.frame(
            x1=xmap[ab[1]], x2=xmap[ab[2]], y=height + i*step, label=sig_pairs$label[k]
          ))
          i <- i + 1L
        }
      }
      
      p <- p + annotate("text", x=0.7, y=ymax_plot*0.98,
                        label=paste0("ANOVA p = ", signif(anov_p, 3)),
                        hjust=0, size=3.6)
    } else {
      note <- paste0(note, " | stats skipped (insufficient n per group)")
    }
  }
  
  if (!is.null(ann) && nrow(ann)) {
    for (k in seq_len(nrow(ann))) {
      p <- p +
        annotate("segment", x=ann$x1[k], xend=ann$x2[k], y=ann$y[k], yend=ann$y[k]) +
        annotate("segment", x=ann$x1[k], xend=ann$x1[k], y=ann$y[k], yend=ann$y[k] - 0.02*ymax_bar) +
        annotate("segment", x=ann$x2[k], xend=ann$x2[k], y=ann$y[k], yend=ann$y[k] - 0.02*ymax_bar) +
        annotate("text", x=(ann$x1[k]+ann$x2[k])/2, y=ann$y[k] + 0.02*ymax_bar,
                 label=ann$label[k], vjust=0)
    }
  }
  
  list(plot=p, note=note, removed=removed)
}

run_pipeline <- function() {
  infile  <- tclvalue(v_infile)
  outbase <- tclvalue(v_outdir)
  outname <- tclvalue(v_outname)
  
  if (!nzchar(infile) || !file.exists(infile)) stop("Input CSV not found.")
  if (!nzchar(outbase) || !dir.exists(outbase)) stop("Output folder (parent) not found.")
  
  group_col <- tclvalue(v_groupcol)
  batch_col <- tclvalue(v_batchcol)
  
  batch_mode <- tclvalue(v_batch_mode)
  group_mode <- tclvalue(v_group_mode)
  
  palette_choice <- tclvalue(v_palette)
  if (!(palette_choice %in% palette_choices)) palette_choice <- "Pastel"
  
  hz <- suppressWarnings(as.numeric(tclvalue(v_hz)))
  if (!is.finite(hz) || hz <= 0) stop("Hz must be a positive number.")
  tol_pct <- suppressWarnings(as.numeric(tclvalue(v_tol_pct)))
  if (!is.finite(tol_pct) || tol_pct < 0 || tol_pct > 100) stop("Tolerance (%) must be between 0 and 100.")
  tol <- tol_pct / 100
  
  expected <- hz * WINDOW_SEC
  lo <- expected * (1 - tol)
  hi <- expected * (1 + tol)
  
  dat <- suppressMessages(readr::read_csv(infile, show_col_types = FALSE))
  if (!(group_col %in% names(dat))) stop("Group column not found: ", group_col)
  if (!("Well" %in% names(dat))) dat$Well <- seq_len(nrow(dat))
  
  has_batch <- batch_col %in% names(dat)
  if (!has_batch) {
    dat$Batch <- factor("Batch1")
    batch_col <- "Batch"
  }
  
  dat <- dat %>%
    rename(Group = all_of(group_col)) %>%
    rename(Batch = all_of(batch_col))
  
  # ---- Critical FIX: force factor types BEFORE droplevels() ----
  dat$Group <- as.factor(dat$Group)
  dat$Batch <- as.factor(dat$Batch)
  
  # Batch selection
  if (batch_mode == "select") {
    keep_batches <- get_listbox_selection(lb_batch, batch_values)
    if (!length(keep_batches)) stop("Batch mode=Select, but no batches selected.")
    dat <- dat %>% filter(as.character(Batch) %in% keep_batches)
    dat$Batch <- droplevels(dat$Batch)
  }
  
  # Group selection + ORDER (uses CURRENT listbox order, after ↑/↓)
  if (group_mode == "select") {
    keep_groups <- get_listbox_selection_in_current_order(lb_group)
    if (!length(keep_groups)) stop("Group mode=Select, but no groups selected.")
    dat <- dat %>% filter(as.character(Group) %in% keep_groups)
    
    # enforce chosen order in plots
    dat$Group <- factor(as.character(dat$Group), levels = keep_groups)
    dat$Group <- droplevels(dat$Group)
  } else {
    # If "All groups", keep dataset ordering (sorted by factor default); no forced preferred order
    dat$Group <- droplevels(dat$Group)
  }
  
  if (!nrow(dat)) stop("No rows left after Batch/Group selection.")
  
  # Peaks QC (always)
  if (!(PEAKS_COL %in% names(dat))) stop("Peaks QC requires column: ", PEAKS_COL)
  dat[[PEAKS_COL]] <- coerce_numeric_like(dat[[PEAKS_COL]])
  
  n_before <- nrow(dat)
  
  removed_peaks <- dat %>%
    filter(is.na(.data[[PEAKS_COL]]) | .data[[PEAKS_COL]] < lo | .data[[PEAKS_COL]] > hi) %>%
    transmute(Well, Group, Batch,
              Num.Peaks = .data[[PEAKS_COL]],
              expected = expected, lo = lo, hi = hi)
  
  dat <- dat %>%
    filter(!is.na(.data[[PEAKS_COL]]),
           .data[[PEAKS_COL]] >= lo,
           .data[[PEAKS_COL]] <= hi)
  
  n_after <- nrow(dat)
  
  safe_message("Peaks QC: Hz=", hz, " window=", WINDOW_SEC, "s → expected=", expected,
               " | tol=", tol_pct, "% → keep [", round(lo,2), ", ", round(hi,2), "]",
               " | removed ", (n_before - n_after), " rows (kept ", n_after, ").")
  
  if (!n_after) stop("All rows removed by peaks QC. Check Hz/tolerance or Num.Peaks.")
  
  dat$Group <- droplevels(dat$Group)
  dat$Batch <- droplevels(dat$Batch)
  
  # Coerce numeric-like metrics
  exclude_cols <- c("Well", PEAKS_COL, "Group", "Batch", "batch", "source_file")
  dat <- as.data.frame(dat)
  to_convert <- setdiff(names(dat), exclude_cols)
  dat[to_convert] <- lapply(dat[to_convert], coerce_numeric_like)
  
  num_cols_all <- setdiff(names(dat)[sapply(dat, is.numeric)], exclude_cols)
  has_data <- vapply(num_cols_all, function(v) sum(is.finite(dat[[v]]), na.rm = TRUE) > 0, logical(1))
  num_cols <- num_cols_all[has_data]
  
  outdir <- file.path(outbase, outname)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  diag_path      <- file.path(outdir, "plot_diagnostics.csv")
  outliers_path  <- file.path(outdir, "outliers_removed.csv")
  peaks_log_path <- file.path(outdir, "peak_filter_log.csv")
  
  readr::write_csv(removed_peaks, peaks_log_path)
  
  if (!length(num_cols)) {
    readr::write_csv(tibble(variable=character(), plotted=logical(), note=character()), diag_path)
    readr::write_csv(tibble(variable=character(), Well=integer(), Group=character(), Batch=character(), value=double()),
                     outliers_path)
    safe_message("No numeric columns with finite values. Logs written to: ", outdir)
    return(invisible(NULL))
  }
  
  diag <- list(); made <- 0L; removed_all <- list()
  
  for (vn in num_cols) {
    res <- try(make_plot_one(dat, vn, palette_choice), silent = TRUE)
    if (is.list(res) && is.null(res$error)) {
      base_name <- file.path(outdir, sanitize_filename(vn))
      save_plot_png_pdf(res$plot, base_name, width = 7, height = 4.6, dpi = 300)
      
      made <- made + 1L
      diag[[length(diag)+1L]] <- data.frame(variable = vn, plotted = TRUE, note = res$note)
      if (!is.null(res$removed) && nrow(res$removed)) removed_all[[length(removed_all)+1L]] <- res$removed
    } else {
      reason <- if (is.list(res)) res$error else as.character(res)
      diag[[length(diag)+1L]] <- data.frame(variable = vn, plotted = FALSE, note = reason)
      if (is.list(res) && !is.null(res$removed) && nrow(res$removed)) removed_all[[length(removed_all)+1L]] <- res$removed
    }
  }
  
  readr::write_csv(bind_rows(diag), diag_path)
  
  if (length(removed_all)) {
    rm_tbl <- bind_rows(removed_all) %>%
      select(variable, Well, Group, Batch, value)
    readr::write_csv(rm_tbl, outliers_path)
  } else {
    readr::write_csv(tibble(variable=character(), Well=integer(), Group=character(), Batch=character(), value=double()),
                     outliers_path)
  }
  
  safe_message("Plots saved (PNG+PDF): ", made, " → ", outdir)
  safe_message("Diagnostics: ", diag_path)
  safe_message("Outliers: ", outliers_path)
  safe_message("Peaks QC log: ", peaks_log_path)
  print(dat %>% count(Group, Batch, name="n"))
  
  invisible(NULL)
}

# RUN button
btn_run <- tkbutton(tt, text="RUN", width=12, command=function(){
  tclvalue(v_status) <- "Running..."
  tryCatch({
    run_pipeline()
    tclvalue(v_status) <- "Done ✓ (check output folder)"
  }, error=function(e){
    tclvalue(v_status) <- paste0("ERROR: ", e$message)
  })
})
tkgrid(btn_run, row=row, column=2, padx=6, pady=8)

tkfocus(tt)
tkwait.window(tt)