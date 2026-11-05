# Acoustic vowel variability: local and global measures
# Run from the repository root: Rscript variability_analysis.R [input.csv]

library(tidyverse)
library(patchwork)
library(ggtext)
library(lmerTest)
library(emmeans)

args <- commandArgs(trailingOnly = TRUE)
input_csv <- if (length(args)) args[[1]] else "data/acoustics_revised_paired.csv"
if (!file.exists(input_csv)) stop("Input CSV not found: ", input_csv)
formant_data <- read_csv(input_csv, show_col_types = FALSE)
required <- c("language", "speaker", "vowel", "word", "token_id", "F1", "F2", "F3")
missing_columns <- setdiff(required, names(formant_data))
if (length(missing_columns)) stop("Missing columns: ", paste(missing_columns, collapse = ", "))
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

set.seed(42)  # Reproducible horizontal jitter

save_figure <- function(plot, name, width, height, units = "cm") {
  ggsave(file.path("outputs/figures", paste0(name, ".png")),
         plot = plot, width = width, height = height,
         units = units, dpi = 400, bg = "white")
}

# 1. Input and speaker normalization ----

formants <- formant_data %>%
  transmute(
    language = language,
    speaker  = speaker,
    vowel    = vowel,
    word     = word,
    token_id = token_id,   # one ID per repetition/token
    F1       = F1,
    F2       = F2,
    F3       = F3
  )

nasal_vowels <- tribble(
  ~language,     ~vowel,
  "French",      "AN",
  "French",      "CN",
  "French",      "EN",
  "French",      "YN",
  "Portuguese",  "EN",
  "Portuguese",  "IN",
  "Portuguese",  "ON",
  "Portuguese",  "UN",
  "Portuguese",  "WN"
)

formants <- formants %>%
  anti_join(
    nasal_vowels,
    by = c("language", "vowel")
  )

formants_analyzed <- formants %>%
  group_by(language, speaker, vowel, word) %>%
  filter(
    sum(complete.cases(F1, F2)) >= 3
  ) %>%
  ungroup() %>%
  filter(
    complete.cases(F1, F2)
  )

write_csv(count(formants_analyzed, language, name = "n_tokens"),
          "outputs/tables/descriptive_token_counts.csv")

# The >=3 complete-token subset above is used for descriptive counts only.
# Preserve the original oral-token input to normalization and metric calculations.
formants <- formants %>%
  group_by(language, speaker) %>%
  mutate(
    F1_z = as.numeric(scale(F1)),
    F2_z = as.numeric(scale(F2)),
    F3_z = as.numeric(scale(F3))
  ) %>%
  ungroup()

# 2. Local acoustic measures ----

#Helper functions
# Mean Euclidean distance from the category centroid
mean_centroid_distance <- function(data, variables) {
  x <- data %>%
    select(all_of(variables)) %>%
    drop_na() %>%
    as.matrix()

  if (nrow(x) < 2) {
    return(NA_real_)
  }

  centroid <- colMeans(x)

  mean(
    sqrt(rowSums(sweep(x, 2, centroid, FUN = "-")^2))
  )
}


# Root-mean-square distance from the centroid
rms_centroid_distance <- function(data, variables) {
  x <- data %>%
    select(all_of(variables)) %>%
    drop_na() %>%
    as.matrix()

  if (nrow(x) < 2) {
    return(NA_real_)
  }

  centroid <- colMeans(x)

  sqrt(
    mean(rowSums(sweep(x, 2, centroid, FUN = "-")^2))
  )
}


# Area of a bivariate p-level data ellipse
ellipse_area <- function(x, y, level = 0.95) {
  xy <- tibble(x = x, y = y) %>%
    drop_na() %>%
    as.matrix()

  if (nrow(xy) < 3) {
    return(NA_real_)
  }

  covariance <- cov(xy)
  determinant <- det(covariance)

  if (!is.finite(determinant) || determinant <= 0) {
    return(NA_real_)
  }

  pi * qchisq(level, df = 2) * sqrt(determinant)
}

##Local acoustic variability by speaker × vowel × word
#This calculates:
#F1, F2, and F3 SD;
#F1, F2, and F3 variance;
#F1, F2, and F3 CV;
#F1–F2 mean centroid distance;
#F1–F2–F3 mean centroid distance;
#RMS centroid distance;
#95% F1–F2 ellipse area.

acoustic_local_word <- formants %>%
  group_by(language, speaker, vowel, word) %>%
  summarise(
    n_tokens = n(),

    # Dimension-specific absolute variability
    F1_mean = mean(F1, na.rm = TRUE),
    F2_mean = mean(F2, na.rm = TRUE),
    F3_mean = mean(F3, na.rm = TRUE),

    F1_sd = sd(F1, na.rm = TRUE),
    F2_sd = sd(F2, na.rm = TRUE),
    F3_sd = sd(F3, na.rm = TRUE),

    F1_variance = var(F1, na.rm = TRUE),
    F2_variance = var(F2, na.rm = TRUE),
    F3_variance = var(F3, na.rm = TRUE),

    # Coefficients of variation
    F1_cv = 100 * F1_sd / F1_mean,
    F2_cv = 100 * F2_sd / F2_mean,
    F3_cv = 100 * F3_sd / F3_mean,

    # Multidimensional variability in normalized space
    mean_distance_F1F2 = mean_centroid_distance(
      pick(everything()),
      c("F1_z", "F2_z")
    ),

    mean_distance_F1F2F3 = mean_centroid_distance(
      pick(everything()),
      c("F1_z", "F2_z", "F3_z")
    ),

    rms_distance_F1F2 = rms_centroid_distance(
      pick(everything()),
      c("F1_z", "F2_z")
    ),

    ellipse_area_F1F2 = ellipse_area(F1_z, F2_z),

    .groups = "drop"
  )


#Average the word-level values to speaker × vowel
acoustic_local_vowel <- acoustic_local_word %>%
  group_by(language, speaker, vowel) %>%
  summarise(
    n_words = n_distinct(word),

    across(
      c(
        F1_sd, F2_sd, F3_sd,
        F1_variance, F2_variance, F3_variance,
        F1_cv, F2_cv, F3_cv,
        mean_distance_F1F2,
        mean_distance_F1F2F3,
        rms_distance_F1F2,
        ellipse_area_F1F2
      ),
      ~ mean(.x, na.rm = TRUE)
    ),

    .groups = "drop"
  ) %>%
  mutate(
    across(
      where(is.numeric),
      ~ ifelse(is.nan(.x), NA_real_, .x)
    )
  )

# 3. Global measures for both English low-vowel choices ----

# Use the same non-English mapping in both versions. The low vowel is
# internally labelled "a" in both cases to keep pairwise metric names stable.
make_iau_centroids <- function(data, english_low) {
  data %>%
    mutate(
      vowel_iau = case_when(
        vowel %in% c("i", "IY") ~ "i",
        vowel %in% c("u", "UW") ~ "u",
        language == "English" & vowel %in% english_low ~ "a",
        language != "English" & vowel %in% c("a", "AA") ~ "a",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(vowel_iau)) %>%
    group_by(language, speaker, vowel_iau) %>%
    summarise(
      F1_centroid = mean(F1_z, na.rm = TRUE),
      F2_centroid = mean(F2_z, na.rm = TRUE),
      .groups = "drop"
    )
}

iau_centroids_AE <- make_iau_centroids(formants, c("AE", "æ"))
iau_centroids_AH <- make_iau_centroids(formants, c("AH", "ɑ"))

#Helper function for triangle area
triangle_area_ <- function(x, y) {
  if (length(x) != 3L || length(y) != 3L) {
    return(NA_real_)
  }

  abs(
    x[1] * (y[2] - y[3]) +
      x[2] * (y[3] - y[1]) +
      x[3] * (y[1] - y[2])
  ) / 2
}

calculate_global_metrics <- function(centroids) {

  # Require three complete, finite vowel centroids per speaker
  complete <- centroids %>%
    filter(
      vowel_iau %in% c("i", "a", "u"),
      is.finite(F1_centroid),
      is.finite(F2_centroid)
    ) %>%
    group_by(language, speaker) %>%
    filter(
      n() == 3L,
      n_distinct(vowel_iau) == 3L
    ) %>%
    ungroup()

  system_metrics <- complete %>%
    group_by(language, speaker) %>%
    summarise(
      triangle_area = triangle_area_(
        F1_centroid, F2_centroid
      ),

      mean_system_radius = mean(
        sqrt(
          (F1_centroid - mean(F1_centroid))^2 +
            (F2_centroid - mean(F2_centroid))^2
        )
      ),

      F1_range = diff(range(F1_centroid)),
      F2_range = diff(range(F2_centroid)),

      .groups = "drop"
    )

  # Self-join the SAME centroid dataset
  pairwise_metrics <- complete %>%
    inner_join(
      complete,
      by = c("language", "speaker"),
      suffix = c("_1", "_2"),
      relationship = "many-to-many"
    ) %>%
    filter(vowel_iau_1 < vowel_iau_2) %>%
    mutate(
      metric = case_when(
        vowel_iau_1 == "a" & vowel_iau_2 == "i" ~ "distance_ia",
        vowel_iau_1 == "a" & vowel_iau_2 == "u" ~ "distance_au",
        vowel_iau_1 == "i" & vowel_iau_2 == "u" ~ "distance_iu"
      ),
      value = sqrt(
        (F1_centroid_1 - F1_centroid_2)^2 +
          (F2_centroid_1 - F2_centroid_2)^2
      )
    ) %>%
    select(language, speaker, metric, value) %>%
    pivot_wider(
      names_from = metric,
      values_from = value
    )

  system_metrics %>%
    left_join(
      pairwise_metrics,
      by = c("language", "speaker")
    )
}

global_metrics_AE <- calculate_global_metrics(iau_centroids_AE)
global_metrics_AH <- calculate_global_metrics(iau_centroids_AH)

# The combined tables contain both system measures and pairwise distances.
global_iau_dispersion_AE <- global_metrics_AE
global_iau_dispersion_AH <- global_metrics_AH

pairwise_long <- function(data) {
  data %>%
    select(language, speaker, `a–i` = distance_ia,
           `a–u` = distance_au, `i–u` = distance_iu) %>%
    pivot_longer(-c(language, speaker), names_to = "contrast", values_to = "distance")
}
iau_pairwise_distances_AE <- pairwise_long(global_metrics_AE)
iau_pairwise_distances_AH <- pairwise_long(global_metrics_AH)

write_csv(acoustic_local_vowel, "outputs/tables/acoustic_local_vowel.csv")
write_csv(global_metrics_AE, "outputs/tables/global_metrics_AE.csv")
write_csv(global_metrics_AH, "outputs/tables/global_metrics_AH.csv")

# 4. Plot settings and local metric figure ----

language_order <- c(
  "Spanish",
  "Portuguese",
  "English",
  "French"
)

language_colors <- c(
  "Spanish" = "#D55E00",
  "Portuguese" = "#009E73",
  "English" = "#0072B2",
  "French" = "#CC79A7"
)

theme_poster <- theme_classic(base_size = 17) +
  theme(
    plot.title = element_text(
      size = 20,
      face = "bold"
    ),
    plot.subtitle = element_text(
      size = 14,
      margin = margin(b = 10)
    ),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    strip.text = element_text(
      size = 15,
      face = "bold"
    ),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13),
    legend.position = "bottom"
  )

#Main local-variability figure
local_metrics <- c(
  "F1_sd",
  "F2_sd",
  "F1_cv",
  "F2_cv",
  "mean_distance_F1F2",
  "ellipse_area_F1F2"
)

metric_labels <- c(
  "F1_sd" = "F1 SD",
  "F2_sd" = "F2 SD",
  "F1_cv" = "F1 CV",
  "F2_cv" = "F2 CV",
  "mean_distance_F1F2" = "Mean F1–F2 distance",
  "ellipse_area_F1F2" = "F1–F2 ellipse area"
)

local_speaker_long <- acoustic_local_vowel %>%
  select(
    language,
    speaker,
    all_of(local_metrics)
  ) %>%
  pivot_longer(
    cols = all_of(local_metrics),
    names_to = "metric",
    values_to = "value"
  ) %>%
  group_by(language, speaker, metric) %>%
  summarise(
    # Equal weighting of the oral vowels within each speaker
    value = mean(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(metric) %>%
  mutate(
    variability_z = as.numeric(
      scale(log(value))
    )
  ) %>%
  ungroup() %>%
  mutate(
    language = factor(
      language,
      levels = language_order
    ),
    metric = factor(
      metric,
      levels = local_metrics,
      labels = metric_labels[local_metrics]
    )
  )

#Calculate descriptive means and 95% normal-approximation intervals
local_summary <- local_speaker_long %>%
  group_by(language, metric) %>%
  summarise(
    mean_z = mean(variability_z),
    se = sd(variability_z) / sqrt(n()),
    lower = mean_z - 1.96 * se,
    upper = mean_z + 1.96 * se,
    .groups = "drop"
  )

#Plot
p_local <- ggplot(
  local_speaker_long,
  aes(
    x = language,
    y = variability_z,
    color = language
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.4,
    color = "grey70"
  ) +
  geom_jitter(
    width = 0.10,
    height = 0,
    alpha = 0.28,
    size = 1.8
  ) +
  geom_errorbar(
    data = local_summary,
    aes(
      x = language,
      ymin = lower,
      ymax = upper,
      color = language
    ),
    inherit.aes = FALSE,
    width = 0.10,
    linewidth = 0.8
  ) +
  geom_point(
    data = local_summary,
    aes(
      x = language,
      y = mean_z,
      color = language
    ),
    inherit.aes = FALSE,
    size = 3.4
  ) +
  facet_wrap(
    ~ metric,
    ncol = 3
  ) +
  scale_color_manual(
    values = language_colors,
    guide = "none"
  ) +
  scale_x_discrete(
    limits = language_order,
    labels = function(x) {
      sprintf(
        "<span style='color: %s;'>%s</span>",
        unname(language_colors[x]),
        x
      )
    }
  ) +
  labs(
    x = NULL,
    y = "Standardized variability"
  ) +
  theme_poster +
  theme(
    axis.text.x = ggtext::element_markdown(
      angle = 30,
      hjust = 1
    )
  )

save_figure(p_local, "local_metric_dependence", 14, 8, "in")

# 5. Global metric figures ----

global_metrics <- c(
  "triangle_area",
  "mean_system_radius",
  "F1_range",
  "F2_range",
  "distance_ia",
  "distance_au",
  "distance_iu"
)

global_metric_labels <- c(
  triangle_area = "Triangle area",
  mean_system_radius = "Mean system radius",
  F1_range = "F1 range",
  F2_range = "F2 range",
  distance_ia = "/i-a/ distance",
  distance_au = "/a-u/ distance",
  distance_iu = "/i-u/ distance"
)

# Standardization is separate within each metric and each vowel-choice version.
plot_global_metrics <- function(global_source, english_low_label) {
  global_speaker_long <- global_source %>%
    select(
      language,
      speaker,
      all_of(global_metrics)
    ) %>%
    pivot_longer(
      cols = all_of(global_metrics),
      names_to = "metric",
      values_to = "value"
    )

  # Log transformation requires finite, positive values.
  # Stop rather than silently remove problematic observations.
  if (any(
    !is.finite(global_speaker_long$value) |
    global_speaker_long$value <= 0
  )) {
    stop("Some global metrics are missing, nonfinite, or nonpositive.")
  }

  global_speaker_long <- global_speaker_long %>%
    group_by(metric) %>%
    mutate(
      dispersion_z = as.numeric(scale(log(value)))
    ) %>%
    ungroup() %>%
    mutate(
      language = factor(language, levels = language_order),
      metric = factor(
        metric,
        levels = global_metrics,
        labels = global_metric_labels[global_metrics]
      )
    )

  # -----------------------------------------------------------
  # 4. Descriptive means and 95% normal-approximation intervals
  # -----------------------------------------------------------

  global_summary <- global_speaker_long %>%
    group_by(language, metric) %>%
    summarise(
      n_speakers = n(),
      mean_z = mean(dispersion_z),
      se = sd(dispersion_z) / sqrt(n()),
      lower = mean_z - 1.96 * se,
      upper = mean_z + 1.96 * se,
      .groups = "drop"
    )

  # -----------------------------------------------------------
  # 5. Plot: same presentation as the local figure
  # -----------------------------------------------------------

  ggplot(
    global_speaker_long,
    aes(
      x = language,
      y = dispersion_z,
      colour = language
    )
  ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.4,
      colour = "grey70"
    ) +
    geom_jitter(
      width = 0.10,
      height = 0,
      alpha = 0.28,
      size = 1.8
    ) +
    geom_errorbar(
      data = global_summary,
      aes(
        x = language,
        ymin = lower,
        ymax = upper,
        colour = language
      ),
      inherit.aes = FALSE,
      width = 0.10,
      linewidth = 0.8
    ) +
    geom_point(
      data = global_summary,
      aes(
        x = language,
        y = mean_z,
        colour = language
      ),
      inherit.aes = FALSE,
      size = 3.4
    ) +
    facet_wrap(~ metric, ncol = 3) +
    scale_colour_manual(
      values = language_colors,
      guide = "none"
    ) +
    labs(
      title = "Cross-linguistic patterns by global variability metric",
      subtitle = paste(
        "English low-vowel representative:", english_low_label
      ),
      x = NULL,
      y = "Standardized variability"
    ) +
    theme_poster +
    theme(
      axis.text.x = element_text(
        angle = 30,
        hjust = 1
      )
    )
}
p_global_AE <- plot_global_metrics(global_metrics_AE, "/æ/")
p_global_AH <- plot_global_metrics(global_metrics_AH, "/ɑ/")
save_figure(p_global_AE, "global_variability_AE", 28, 24)
save_figure(p_global_AH, "global_variability_AH", 28, 24)

# 6. English low-vowel comparison ----

#Global triangle area with both English choices
global_triangle_plot <- bind_rows(
  global_iau_dispersion_AE %>%
    filter(language != "English") %>%
    mutate(system = language),

  global_iau_dispersion_AE %>%
    filter(language == "English") %>%
    mutate(system = "English: /i æ u/"),

  global_iau_dispersion_AH %>%
    filter(language == "English") %>%
    mutate(system = "English: /i ɑ u/")
) %>%
  mutate(
    system = factor(
      system,
      levels = c(
        "Spanish",
        "Portuguese",
        "English: /i æ u/",
        "English: /i ɑ u/",
        "French"
      )
    )
  )

#summarize
global_triangle_summary <- global_triangle_plot %>%
  group_by(system) %>%
  summarise(
    mean = mean(triangle_area),
    se = sd(triangle_area) / sqrt(n()),
    lower = mean - 1.96 * se,
    upper = mean + 1.96 * se,
    .groups = "drop"
  )

#plot
p_triangle <- ggplot(
  global_triangle_plot,
  aes(
    x = system,
    y = triangle_area
  )
) +
  geom_jitter(
    width = 0.10,
    height = 0,
    alpha = 0.30,
    size = 2,
    color = "grey35"
  ) +
  geom_errorbar(
    data = global_triangle_summary,
    aes(
      x = system,
      ymin = lower,
      ymax = upper
    ),
    inherit.aes = FALSE,
    width = 0.10,
    linewidth = 0.9,
    color = "#0072B2"
  ) +
  geom_point(
    data = global_triangle_summary,
    aes(
      x = system,
      y = mean
    ),
    inherit.aes = FALSE,
    size = 4,
    color = "#0072B2"
  ) +
  labs(
    title = "English global dispersion depends on the selected low vowel",
    subtitle = "Points show speakers; large points and error bars show means and 95% intervals",
    x = NULL,
    y = "Triangle area"
  ) +
  theme_poster +
  theme(
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    )
  )

p_triangle

#Speaker-by-speaker AE–AH comparison
#This figure shows that the English change is systematic
#across speakers and that its direction differs by measure.
english_global_comparison <- global_iau_dispersion_AE %>%
  filter(language == "English") %>%
  select(
    speaker,
    triangle_area,
    mean_system_radius,
    F1_range,
    F2_range
  ) %>%
  pivot_longer(
    cols = -speaker,
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(low_vowel = "AE: /æ/") %>%
  bind_rows(
    global_iau_dispersion_AH %>%
      filter(language == "English") %>%
      select(
        speaker,
        triangle_area,
        mean_system_radius,
        F1_range,
        F2_range
      ) %>%
      pivot_longer(
        cols = -speaker,
        names_to = "metric",
        values_to = "value"
      ) %>%
      mutate(low_vowel = "AH: /ɑ/")
  ) %>%
  mutate(
    low_vowel = factor(
      low_vowel,
      levels = c("AE: /æ/", "AH: /ɑ/")
    ),
    metric = factor(
      metric,
      levels = c(
        "triangle_area",
        "mean_system_radius",
        "F1_range",
        "F2_range"
      ),
      labels = c(
        "Triangle area",
        "Mean radius",
        "F1 range",
        "F2 range"
      )
    )
  )

#plot
p_english_paired <- ggplot(
  english_global_comparison,
  aes(
    x = low_vowel,
    y = value,
    group = speaker
  )
) +
  geom_line(
    color = "grey65",
    linewidth = 0.7,
    alpha = 0.8
  ) +
  geom_point(
    aes(color = low_vowel),
    size = 2.5
  ) +
  facet_wrap(
    ~ metric,
    scales = "free_y",
    ncol = 2
  ) +
  scale_color_manual(
    values = c(
      "AE: /æ/" = "#E69F00",
      "AH: /ɑ/" = "#0072B2"
    ),
    guide = "none"
  ) +
  labs(
    title = "Changing the English low vowel alters global dispersion",
    subtitle = "Lines connect the two estimates for each speaker",
    x = "English low vowel",
    y = "Dispersion"
  ) +
  theme_poster

p_english_paired

#Which triangle side produces the change?
english_pairwise <- iau_pairwise_distances_AE %>%
  filter(language == "English") %>%
  mutate(low_vowel = "AE: /æ/") %>%
  bind_rows(
    iau_pairwise_distances_AH %>%
      filter(language == "English") %>%
      mutate(low_vowel = "AH: /ɑ/")
  ) %>%
  mutate(
    low_vowel = factor(
      low_vowel,
      levels = c("AE: /æ/", "AH: /ɑ/")
    ),
    contrast = recode(
      contrast,
      "a–i" = "low–i",
      "a–u" = "low–u",
      "i–u" = "i–u"
    )
  )

#plot
p_pairwise <- ggplot(
  english_pairwise,
  aes(
    x = low_vowel,
    y = distance,
    group = speaker
  )
) +
  geom_line(
    color = "grey70",
    linewidth = 0.6
  ) +
  geom_point(
    aes(color = low_vowel),
    size = 2.3
  ) +
  facet_wrap(
    ~ contrast,
    nrow = 1
  ) +
  scale_color_manual(
    values = c(
      "AE: /æ/" = "#E69F00",
      "AH: /ɑ/" = "#0072B2"
    ),
    guide = "none"
  ) +
  labs(
    title = "The English difference is concentrated in the low–u distance",
    x = "English low vowel",
    y = "Pairwise distance"
  ) +
  theme_poster

p_pairwise

#combined figure
p_global_combined <-
  p_triangle /
  (p_english_paired + p_pairwise) +
  plot_annotation(
    tag_levels = "A"
  )

save_figure(p_triangle, "triangle_area_both_english_vowels", 38, 14)
save_figure(p_global_combined, "english_global_dispersion", 17, 13, "in")

# Overlaid vowel-space triangles (equal weighting of words, then speakers).
en <- formants %>%
  filter(
    language == "English",
    vowel %in% c("IY", "AE", "AH", "UW")
  ) %>%
  transmute(
    speaker,
    vowel = as.character(vowel),
    word,
    F1 = F1_z,
    F2 = F2_z
  ) %>%
  filter(is.finite(F1), is.finite(F2))

# Equal weighting of words within speaker and vowel
speaker_centroids <- en %>%
  group_by(speaker, vowel, word) %>%
  summarise(
    F1 = mean(F1),
    F2 = mean(F2),
    .groups = "drop"
  ) %>%
  group_by(speaker, vowel) %>%
  summarise(
    F1 = mean(F1),
    F2 = mean(F2),
    .groups = "drop"
  )

# Retain the same speakers for both triangles
complete_speakers <- speaker_centroids %>%
  count(speaker) %>%
  filter(n == 4)

speaker_centroids <- speaker_centroids %>%
  semi_join(complete_speakers, by = "speaker")

# Equal weighting of speakers
centroids <- speaker_centroids %>%
  group_by(vowel) %>%
  summarise(
    F1 = mean(F1),
    F2 = mean(F2),
    .groups = "drop"
  )

make_triangle <- function(low, label) {
  centroids %>%
    filter(vowel %in% c("IY", low, "UW")) %>%
    mutate(
      order = match(vowel, c("IY", low, "UW")),
      triangle = label
    ) %>%
    arrange(order)
}

triangles <- bind_rows(
  make_triangle("AE", "/i–æ–u/"),
  make_triangle("AH", "/i–ɑ–u/")
) %>%
  mutate(triangle = factor(
    triangle, levels = c("/i–æ–u/", "/i–ɑ–u/")
  ))

# Put range markers beneath the vowel space.
# With a reversed F1 axis, larger y values appear lower.
span <- max(diff(range(centroids$F1)), 0.5)
bottom <- max(centroids$F1)

ranges <- triangles %>%
  group_by(triangle) %>%
  summarise(
    xmin = min(F2),
    xmax = max(F2),
    .groups = "drop"
  ) %>%
  mutate(
    y = bottom + c(0.28, 0.52) * span,
    label = paste0(triangle, " F2 range")
  )

ipa <- c(IY = "i", AE = "æ", AH = "ɑ", UW = "u")
colours <- c("/i–æ–u/" = "#0072B2",
             "/i–ɑ–u/" = "#D55E00")

plot_triangles_EN <- ggplot() +
  geom_polygon(
    data = triangles,
    aes(F2, F1, group = triangle,
        colour = triangle, fill = triangle),
    alpha = 0.10,
    linewidth = 1
  ) +
  geom_point(
    data = centroids,
    aes(F2, F1),
    size = 2,
    colour = "grey25"
  ) +
  geom_text(
    data = centroids,
    aes(F2, F1, label = ipa[vowel]),
    nudge_y = -0.065 * span,
    size = 6
  ) +
  geom_segment(
    data = ranges,
    aes(x = xmin, xend = xmax,
        y = y, yend = y, colour = triangle),
    linewidth = 1.2
  ) +
  geom_segment(
    data = ranges,
    aes(x = xmin, xend = xmin,
        y = y - 0.035 * span,
        yend = y + 0.035 * span,
        colour = triangle)
  ) +
  geom_segment(
    data = ranges,
    aes(x = xmax, xend = xmax,
        y = y - 0.035 * span,
        yend = y + 0.035 * span,
        colour = triangle)
  ) +
  geom_text(
    data = ranges,
    aes(x = (xmin + xmax) / 2,
        y = y + 0.085 * span,
        label = label, colour = triangle),
    size = 4
  ) +
  scale_x_reverse(expand = expansion(mult = 0.12)) +
  scale_y_reverse(expand = expansion(mult = 0.12)) +
  scale_colour_manual(values = colours, name = NULL) +
  scale_fill_manual(values = colours, name = NULL) +
  coord_equal() +
  labs(
    x = "F2 (speaker-standardized)",
    y = "F1 (speaker-standardized)"
  ) +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )

save_figure(plot_triangles_EN, "english_triangles", 22, 17)

# 7. Mixed-effects models ----

local_model_data <- acoustic_local_vowel %>%
  select(
    language,
    speaker,
    vowel,
    all_of(local_metrics)
  ) %>%
  pivot_longer(
    cols = all_of(local_metrics),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(
    is.finite(value),
    value > 0
  ) %>%
  group_by(metric) %>%
  mutate(
    value_log = log(value),
    variability_z = as.numeric(scale(value_log))
  ) %>%
  ungroup() %>%
  mutate(
    language = factor(
      language,
      levels = c(
        "Spanish",
        "Portuguese",
        "English",
        "French"
      )
    ),
    metric = factor(
      metric,
      levels = local_metrics
    ),
    speaker_vowel = interaction(
      language,
      speaker,
      vowel,
      drop = TRUE
    )
  )

model_metric_dependence <- lmer(
  variability_z ~ language * metric +
    (1 | speaker) +
    (1 | speaker_vowel),
  data = local_model_data,
  REML = FALSE
)

local_anova <- anova(model_metric_dependence, type = 3)
print(local_anova)

local_emmeans <- emmeans(
  model_metric_dependence,
  ~ language | metric
)

local_language_pairs <- pairs(
  local_emmeans,
  adjust = "tukey"
)


write_csv(as.data.frame(summary(local_emmeans)), "outputs/tables/local_emmeans.csv")
write_csv(as.data.frame(summary(local_language_pairs)), "outputs/tables/local_language_pairs.csv")

# Preserve the AE choice from the supplied global model.
global_source <- global_metrics_AE


global_model_data <- global_source %>%
  select(
    language,
    speaker,
    all_of(global_metrics)
  ) %>%
  pivot_longer(
    cols = all_of(global_metrics),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(
    is.finite(value),
    value > 0
  ) %>%
  group_by(metric) %>%
  mutate(
    value_log = log(value),
    dispersion_z = as.numeric(scale(value_log))
  ) %>%
  ungroup() %>%
  mutate(
    language = factor(
      language,
      levels = c(
        "Spanish",
        "Portuguese",
        "English",
        "French"
      )
    ),
    metric = factor(
      metric,
      levels = global_metrics
    ),
    speaker_id = interaction(
      language,
      speaker,
      drop = TRUE
    )
  )

# Sum contrasts for Type III tests
contrasts(global_model_data$language) <- contr.sum(
  nlevels(global_model_data$language)
)

contrasts(global_model_data$metric) <- contr.sum(
  nlevels(global_model_data$metric)
)

# Fit model
model_global_metric_dependence <- lmer(
  dispersion_z ~ language * metric +
    (1 | speaker_id),
  data = global_model_data,
  REML = FALSE
)

# Type III tests with Satterthwaite degrees of freedom
global_anova <- anova(
  model_global_metric_dependence,
  type = 3,
  ddf = "Satterthwaite"
)

# Language comparisons within each metric
global_emmeans <- emmeans(
  model_global_metric_dependence,
  ~ language | metric,
  lmer.df = "satterthwaite"
)

global_language_pairs <- pairs(global_emmeans, adjust = "tukey")
write_csv(as.data.frame(summary(global_emmeans)), "outputs/tables/global_emmeans_AE.csv")
write_csv(as.data.frame(summary(global_language_pairs)), "outputs/tables/global_language_pairs_AE.csv")

print(global_anova)

save_model_checks <- function(model, anova_table, label) {
  write_csv(tibble::rownames_to_column(as.data.frame(anova_table), "term"),
            file.path("outputs/tables", paste0(label, "_anova.csv")))
  grDevices::pdf(file.path("outputs/figures", paste0(label, "_diagnostics.pdf")))
  on.exit(grDevices::dev.off())
  print(plot(model))
  qqnorm(resid(model))
  qqline(resid(model))
  message(label, ": singular fit = ", lme4::isSingular(model, tol = 1e-4))
}

save_model_checks(model_metric_dependence, local_anova, "local")
save_model_checks(model_global_metric_dependence, global_anova, "global_AE")

# 8. Methods schematic using simulated coordinates ----

set.seed(42)

blue   <- "#0072B2"
orange <- "#D55E00"
green  <- "#009E73"

# -----------------------------------------------------------
# 1. LOCAL VARIABILITY: simulated tokens from one category
# -----------------------------------------------------------

n <- 45
a <- rnorm(n)
b <- rnorm(n)

tokens <- tibble(
  F2 = 0.55 * a,
  F1 = 0.28 * (0.35 * a + sqrt(1 - 0.35^2) * b)
)

centre <- tokens %>%
  summarise(F1 = mean(F1), F2 = mean(F2))

mu <- c(centre$F2, centre$F1)
S  <- cov(tokens[c("F2", "F1")])

# Explicit 95% Gaussian distribution ellipse
theta <- seq(0, 2 * pi, length.out = 200)
circle <- cbind(cos(theta), sin(theta))

ellipse_xy <- sqrt(qchisq(0.95, df = 2)) *
  circle %*% chol(S)

ellipse_xy <- sweep(ellipse_xy, 2, mu, "+")

ellipse <- tibble(
  F2 = ellipse_xy[, 1],
  F1 = ellipse_xy[, 2]
)

# Illustrate distances for a few tokens only
distance_lines <- tokens %>%
  slice(c(3, 9, 17, 26, 35)) %>%
  mutate(F2_end = centre$F2, F1_end = centre$F1)

sd_f1 <- sd(tokens$F1)
sd_f2 <- sd(tokens$F2)

# SD markers positioned outside the token cloud
x_sd <- 1.65
y_sd <- 1.10
cap  <- 0.045

p_methods_local <- ggplot() +
  geom_polygon(
    data = ellipse, aes(F2, F1),
    fill = blue, alpha = 0.12,
    colour = blue, linewidth = 0.7
  ) +
  geom_segment(
    data = distance_lines,
    aes(x = F2, y = F1, xend = F2_end, yend = F1_end),
    colour = orange, linewidth = 0.65
  ) +
  geom_point(
    data = tokens, aes(F2, F1),
    colour = blue, alpha = 0.7, size = 1.8
  ) +
  geom_point(
    data = centre, aes(F2, F1),
    shape = 4, stroke = 1.3, size = 7
  ) +

  # Horizontal marker: F2 mean ±1 SD
  annotate(
    "segment",
    x = mu[1] - sd_f2, xend = mu[1] + sd_f2,
    y = y_sd, yend = y_sd,
    linewidth = 0.8
  ) +
  annotate(
    "segment",
    x = c(mu[1] - sd_f2, mu[1] + sd_f2),
    xend = c(mu[1] - sd_f2, mu[1] + sd_f2),
    y = y_sd - cap, yend = y_sd + cap
  ) +
  annotate(
    "text", x = mu[1], y = y_sd + 0.22,
    label = "F2: mean ±1 SD", size = 6
  ) +

  # Vertical marker: F1 mean ±1 SD
  annotate(
    "segment",
    x = x_sd, xend = x_sd,
    y = mu[2] - sd_f1, yend = mu[2] + sd_f1,
    linewidth = 0.8
  ) +
  annotate(
    "segment",
    x = x_sd - cap, xend = x_sd + cap,
    y = c(mu[2] - sd_f1, mu[2] + sd_f1),
    yend = c(mu[2] - sd_f1, mu[2] + sd_f1)
  ) +
  annotate(
    "text", x = x_sd + 0.25, y = mu[2],
    label = "F1: mean ±1 SD",
    angle = 90, size = 6
  ) +
  annotate(
    "text",
    x = -0.65, y = -0.45,
    label = "95% ellipse area",
    colour = blue, size = 7
  ) +

  labs(
    title = "Local variability",
    subtitle = "Variation within a vowel category",
    caption = paste(
      "Orange lines: token-to-centroid distances",
      "× Centroid",
      sep = "\n"
    ),
    x = "F2", y = "F1"
  )

# -----------------------------------------------------------
# 2. GLOBAL DISPERSION: three simulated category centroids
# -----------------------------------------------------------

vowels <- tibble(
  vowel = c("i", "a", "u"),
  F2 = c(1.25, 0.05, -1.10),
  F1 = c(-0.80, 0.90, -0.65)
)

global_centre <- vowels %>%
  summarise(F2 = mean(F2), F1 = mean(F1))

radius_lines <- vowels %>%
  mutate(
    F2_end = global_centre$F2,
    F1_end = global_centre$F1
  )

range_min <- min(vowels$F2)
range_max <- max(vowels$F2)
range_y   <- 1.30

p_methods_global <- ggplot() +
  geom_polygon(
    data = vowels, aes(F2, F1),
    fill = blue, alpha = 0.12,
    colour = blue, linewidth = 0.8
  ) +
  geom_segment(
    data = radius_lines,
    aes(x = F2, y = F1, xend = F2_end, yend = F1_end),
    colour = green, linetype = "dashed", linewidth = 0.8
  ) +
  geom_point(
    data = vowels, aes(F2, F1),
    size = 2.5, colour = blue
  ) +
  geom_text(
    data = vowels, aes(F2, F1, label = vowel),
    nudge_x = 0.13, nudge_y = -0.10, size = 8
  ) +
  geom_point(
    data = global_centre, aes(F2, F1),
    shape = 4, stroke = 1.3, size = 7
  ) +
  annotate(
    "text", x = 0.05, y = 0.45,
    label = "Triangle area", size = 7, colour = blue
  ) +

  # F2 range bracket
  annotate(
    "segment",
    x = range_min, xend = range_max,
    y = range_y, yend = range_y,
    linewidth = 0.8
  ) +
  annotate(
    "segment",
    x = c(range_min, range_max),
    xend = c(range_min, range_max),
    y = range_y - cap, yend = range_y + cap
  ) +
  annotate(
    "text",
    x = (range_min + range_max) / 2,
    y = range_y + 0.20,
    label = "F2 range", size = 6
  ) +
  labs(
    title = "Global dispersion",
    subtitle = "Configuration of vowel-category centroids",
    caption = paste(
      "Green lines: vowel-to-centre distances",
      "× System centre",
      sep = "\n"
    ),
    x = "F2", y = "F1"
  )

# -----------------------------------------------------------
# 3. Matching scales and combined figure
# -----------------------------------------------------------

style_panel <- function(p) {
  p +
    scale_x_reverse(limits = c(2.15, -1.75)) +
    scale_y_reverse(limits = c(1.75, -1.45)) +
    coord_equal() +
    theme_classic(base_size = 17) +
    theme(
      plot.title = element_text(face = "bold", size = 20),
      plot.subtitle = element_text(size = 15),
      plot.caption = element_text(
        hjust = 0, size = 14, lineheight = 1.2
      ),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = margin(10, 12, 10, 12)
    )
}

plot_methods <- (
  style_panel(p_methods_local) | style_panel(p_methods_global)
) +
  plot_annotation(
    caption = "Schematic illustration using simulated coordinates"
  )

save_figure(plot_methods, "variability_methods", 28, 15)
writeLines(capture.output(sessionInfo()), "outputs/session-info.txt")
