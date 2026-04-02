<div align="center">

# 🧬 qPCR Analysis Pipeline
### Host-DNA Depletion & Enrichment Broth Optimisation

*R Markdown workflows for foodborne pathogen detection method development*

---

![R](https://img.shields.io/badge/R-%E2%89%A54.2.0-276DC3?style=flat-square&logo=r&logoColor=white)
![License](https://img.shields.io/badge/License-Academic-8B5CF6?style=flat-square)
![Status](https://img.shields.io/badge/Status-PhD%20Thesis%20Chapter%203-10B981?style=flat-square)
![Targets](https://img.shields.io/badge/Targets-Campylobacter%20%7C%20Salmonella-F59E0B?style=flat-square)

</div>

---

## 📌 Overview

This repository contains six R Markdown analysis scripts developed for **Chapter 3 (Method Development)** of a PhD thesis. The scripts evaluate two key methodological questions in the detection of foodborne pathogens from chicken samples:

| Aim | Description |
|-----|-------------|
| 🔬 **Host-DNA Depletion** | Does a depletion protocol reduce chicken (host) DNA while retaining *Campylobacter* and *Salmonella* target DNA across multiple sample preparation methods? |
| 🧫 **Enrichment Optimisation** | Which selective enrichment broth best supports growth of *Campylobacter* and *Salmonella* over 24 hours? |

---

## 📂 Repository Structure

```
📦 qPCR-method-development/
│
├── 🐓 Host-DNA Depletion
│   ├── qPCR_chicken_host-DNA_depletion_enriched.Rmd
│   ├── qPCR_chicken_host-DNA_depletion_4cheesecloth_samples.rmd
│   ├── qPCR_Campylobacter_host-DNA_dep.rmd
│   ├── qPCR_Salmonella_host-DNA_depletion.rmd
│   └── qPCR_filtered_analysis.rmd
│
└── 🧫 Enrichment Optimisation
    └── qPCR_graphs_enrichment.Rmd
```

---

## 📋 Script Descriptions

### 🐓 Host-DNA Depletion Scripts

<details>
<summary><strong>qPCR_chicken_host-DNA_depletion_enriched.Rmd</strong> — Chicken DNA in enriched & filtered samples</summary>

Analyses chicken (host) DNA before and after depletion across cheesecloth-filtered, centrifugation + membrane-filtered, filter-only, RVS-enriched, and CAT-enriched samples.

- Paired t-tests and Wilcoxon signed-rank tests
- Assumption checking (Shapiro-Wilk; IQR outlier rule)
- Cohen's d effect sizes
- Exports a complete PDF report and Excel summary

</details>

<details>
<summary><strong>qPCR_chicken_host-DNA_depletion_4cheesecloth_samples.rmd</strong> — Cheesecloth & filtered samples with non-detect handling</summary>

Focuses on cheesecloth and membrane-filtered samples with explicit handling of non-detects via three complementary statistical approaches.

- Approach 1: Exclusion of non-detect pairs (conservative)
- Approach 2: Detection-limit substitution (sensitivity analysis)
- Approach 3: Sign test (direction-of-change for non-detects)
- Combined Ct + concentration panel plots with non-detect flagging (×)

</details>

<details>
<summary><strong>qPCR_Campylobacter_host-DNA_dep.rmd</strong> — <em>Campylobacter</em> DNA across non-enriched, RVS & CAT samples</summary>

Assesses whether host-DNA depletion affects *Campylobacter* DNA recovery across three sample types.

- Paired tests selected automatically by data structure
- Three statistical approaches for non-detects
- Faceted panel plots per sample type
- Final summary table of depletion statistics

</details>

<details>
<summary><strong>qPCR_Salmonella_host-DNA_depletion.rmd</strong> — <em>Salmonella</em> DNA across non-enriched, RVS & CAT samples</summary>

Mirrors the *Campylobacter* script for *Salmonella* across the same sample types.

- Sign test, paired t-test, and detection-limit substitution
- Exports results to Excel
- Consistent colour schemes across all Salmonella plots

</details>

<details>
<summary><strong>qPCR_filtered_analysis.rmd</strong> — All three DNA targets in filtered samples</summary>

Analyses chicken host DNA, *Campylobacter* DNA, and *Salmonella* DNA within a single unified script for centrifugation + filter and filter-only samples.

- Shared statistical function applied across all three organisms
- Separate high-resolution panel plots per organism (600 dpi)
- Automatic test selection (Wilcoxon vs paired t-test)

</details>

---

### 🧫 Enrichment Optimisation Script

<details>
<summary><strong>qPCR_graphs_enrichment.Rmd</strong> — Broth comparison for <em>Campylobacter</em> & <em>Salmonella</em></summary>

Compares four enrichment broths over a 0–24 hour time course for two pathogen targets.

**Broths tested:**

| Broth | Abbreviation |
|-------|-------------|
| BF Bolton Base | BF-BB |
| BF Bolton Base + CAT supplement | CAT |
| Buffered Peptone Water | BPW |
| Rappaport-Vassiliadis Soya | RVS |

**Outputs:**
- Time-course line plots (mean ± SE with individual replicates)
- Baseline vs 24-hour boxplots per broth
- Between-broth comparison plots at 24 hours
- Dual-axis plots (Ct value + DNA concentration)
- Wilcoxon rank-sum and Kruskal-Wallis tests
- Fold change calculated from both ΔCt (2^-ΔCt) and absolute concentration

</details>

---

## 🧪 Methods at a Glance

### Sample Types

| Category | Sample Types |
|----------|-------------|
| Pre-enrichment | Non-enriched, Cheesecloth-filtered |
| Filtered | Centrifugation + membrane filter, Membrane filter only |
| Enriched | RVS broth, CAT (Bolton Broth + CAT supplement) |

### qPCR Targets

| Target | Role |
|--------|------|
| Chicken 18S rRNA gene | Host DNA (to be depleted) |
| *Campylobacter* spp. | Pathogen target (to be retained) |
| *Salmonella* spp. | Pathogen target (to be retained) |

### Statistical Framework

```
Host-DNA Depletion                    Enrichment Optimisation
─────────────────────────             ─────────────────────────────
Primary: Paired t-test (1-tailed)     Within-broth: Wilcoxon rank-sum
Alternative: Wilcoxon signed-rank     Between-broth: Kruskal-Wallis
Supplementary: Sign test              Fold change: 2^(-ΔCt) & [conc₂₄/conc₀]
Effect size: Cohen's d, fold change   Significance threshold: α = 0.05
Non-detects: 3 complementary methods
```

---

## ⚙️ Dependencies

```r
install.packages(c(
  "tidyverse",   # Data wrangling & plotting
  "ggplot2",     # Visualisation
  "ggpubr",      # Publication-ready plots
  "rstatix",     # Statistical tests
  "readxl",      # Import .xlsx data
  "patchwork",   # Multi-panel figures
  "gridExtra",   # Grid layouts
  "broom",       # Tidy model outputs
  "knitr",       # Report generation
  "kableExtra",  # Formatted tables
  "scales",      # Axis scaling
  "writexl"      # Export to Excel
))
```

> R version ≥ 4.2.0 recommended.

---

## 📊 Outputs

Each script generates:

- 📈 **High-resolution PNG plots** (600 dpi) — all figures
- 📋 **Summary statistics tables** — rendered in HTML or PDF
- 📁 **Excel exports** — statistical results (selected scripts)
- 📄 **PDF reports** — complete analysis summary (selected scripts)

---

## 💾 Data

Raw data are stored as `.xlsx` files and referenced by absolute path within each script. To reproduce the analyses, update the `read_excel()` paths to point to your local copies of the data files.

---

## 📖 Reference

These scripts were developed as part of an unpublished PhD thesis:

> [Author Surname, Initials.] (*in preparation*). *[Thesis Title]*. PhD thesis, [University Name].

*This repository supports the analyses reported in Chapter 3 (Method Development).*

---

## 📄 License

This code is shared for academic reproducibility purposes. Please contact the author before reuse in published work.

---

<div align="center">

*Developed for PhD thesis method development — Chapter 3*

</div>
