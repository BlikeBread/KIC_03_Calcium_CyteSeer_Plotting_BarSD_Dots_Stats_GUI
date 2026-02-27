# KIC Spheroid Calcium Pipeline (Script 03/03)

This repository contains an R script for generating **publication-ready bar plots with statistical analysis** from merged calcium datasets generated in **Script 02a or Script 02b** of the KIC Spheroid Calcium pipeline.

This script represents the third and final step (**03/03**) of a modular **KIC Spheroid Calcium analysis pipeline**.

The pipeline applies strict peak-quality filtering, performs standardized outlier detection, and generates reproducible statistical visualizations exported in both PNG and PDF formats.

---

## What the pipeline does

Starting from a **merged dataset (`merged.csv` or `merged_raw.csv`)** generated in Script 02a or 02b, respectively, the script:

- Prompts the user to select input file and output folder via GUI (`tcltk`)  
- Allows optional selection of:
  - Specific `Batch` values  
  - Specific experimental `Group` values  
- Allows user-defined **group ordering** (↑ / ↓ reordering in GUI)  
- Applies strict peak-quality control:
  - Uses fixed `Num.Peaks` column  
  - Validates expected peaks = Hz × 10-second window  
  - Applies user-defined tolerance (%)  
- Performs automatic numeric coercion and column harmonization  
- Applies within-group outlier detection (1.5×IQR rule)  
- Generates standardized bar plots:
  - Mean ± SD bars  
  - Individual well dots  
  - Batch-based color coding  
- Performs statistical testing:
  - Welch t-test (2 groups)  
  - One-way ANOVA (≥3 groups)  
  - Tukey post-hoc comparisons  
  - Automatic significance annotation (*, **, ***, ****)  
- Exports each metric as:
  - High-resolution PNG  
  - Vector PDF  
- Generates diagnostic and filtering reports  

All file paths are selected interactively via GUI dialogs (`tcltk`).

---

## Required inputs

### Step 02a or 02b output

- `merged.csv` (single-batch, Script 02a)  
**or**
- `merged_raw.csv` (multi-batch, Script 02b)

Required columns:

- `Well`  
- `Group`  
- `Num.Peaks`  

If multi-batch:

- `Batch`  

Optional but recommended:

- `RowGroup`  
- `keep_1hz`  
- `valid_well`  

---

## Generated outputs

The script generates an output folder containing:

- `*.png` — High-resolution statistical plots  
- `*.pdf` — Vector publication-ready plots  
- `plot_diagnostics.csv` — Per-metric plotting summary  
- `outliers_removed.csv` — IQR-filtered values  
- `peak_filter_log.csv` — Peak filtering report  

Each plot is standardized in layout and formatting to ensure cross-experiment consistency.

---

## Typical use cases

- Final statistical comparison across experimental groups  
- Generation of publication-ready figures  
- Batch-aware calcium signal visualization  
- Quality-controlled statistical reporting  
- Reproducible data visualization for manuscripts and reports  

---

## Position in the KIC Pipeline

This script is **Script 03** of a structured workflow:

- Script 01 – Raw CyteSeer CSV processing and MEANS summary generation  
- Script 02a – Single-batch MEANS merging and filtering  
- Script 02b – Cross-batch merging  
- **Script 03 (this repository)** – Statistical analysis and visualization  

---

## Methods Description

Merged calcium transient MEANS datasets were processed using a custom R-based GUI pipeline. Peak-quality control was applied based on expected frequency thresholds, and within-group outlier detection was performed using the 1.5×IQR rule. Statistical comparisons were conducted using Welch t-tests or one-way ANOVA with Tukey post-hoc correction where appropriate. Publication-ready figures were generated and exported in both raster and vector formats to ensure reproducibility and consistency across experiments.

---

## Authorship

This script was developed by **Michele Buono, Talitha Spanjersberg, Nikki Scheen, Nina van der Wilt** and can be used freely for research purposes, provided appropriate citation of the authors.

**Buono, M. F., Spanjersberg, T., Scheen, N., & van der Wilt, C. N. (2026). KIC Spheroid Calcium Pipeline (Script 02b/03) – Multi-Batch Dataset Integration (v.1.0). Zenodo. https://doi.org/10.5281/zenodo.18799739**

The overall workflow, structure, and clarity of the pipeline were iteratively refined with assistance from **OpenAI – ChatGPT 5.2**, which was used as a tool to improve code organization, documentation, and usability.
