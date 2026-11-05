# Acoustic vowel variability across four languages

Supplementary R code for the poster analysis of local acoustic variability and global vowel-space dispersion in Spanish, Portuguese, English, and French.

## Run the analysis

Install the R packages once:

```r
install.packages(c("tidyverse", "patchwork", "ggtext", "lme4", "lmerTest", "emmeans"))
```

From the repository root, put the input at `data/acoustics_revised_paired.csv` and run:

```sh
Rscript variability_analysis.R
```

Alternatively, pass a path to the input CSV as the first argument. The input is not bundled. The script expects one row per token with columns `language`, `speaker`, `vowel`, `word`, `token_id`, `F1`, `F2`, and `F3`. Formants are in Hz; language names are English, French, Portuguese, and Spanish. The code assumes the source dataset's token selection and vowel coding. A separate count of groups with at least three complete F1/F2 observations is descriptive only and does not filter the metric inputs.

## Analysis and outputs

- Local measures are calculated within speaker, vowel, and word, then averaged across words within speaker and vowel. SD and CV use Hz; multidimensional measures use speaker-standardized formants. The local figure averages across vowels within each speaker; the local model retains speaker–vowel observations.
- Global metrics use speaker-specific centroids from normalized formants. The English low vowel is /æ/ (AE) or /ɑ/ (AH); the other languages use the same point-vowel mapping in both versions. The global model uses AE, as specified in the supplied script.
- Plotted cross-language values are log-transformed and standardized separately by metric across the pooled languages. Each global figure has its own reference distribution. Use the original global metric values for direct AE–AH comparisons.
- Results are written to `outputs/tables/`; plots and model diagnostics to `outputs/figures/`. Package versions for each run are recorded in `outputs/session-info.txt`.

| Figure | Filename |
|---|---|
| Local metric comparison | `figures/local_metric_dependence.png` |
| Global metric comparisons | `global_variability_AE.png`, `global_variability_AH.png` |
| Triangle area and paired English comparisons | `english_global_dispersion.png` |
| English vowel-space triangle overlay | `english_triangles.png` |
| Methods schematic (simulated coordinates) | `variability_methods.png` |

)
# vowel-variability-supplement
R code and supplementary figures for a cross-linguistic study of acoustic vowel variability.

## Supplementary figures

- [Global dispersion metrics with English /æ/](figures/global_variability_AE.png)
- [Global dispersion metrics with English /ɑ/](figures/global_variability_AH.png)
- [Triangle area and paired English comparisons](figures/english_global_dispersion.png)
