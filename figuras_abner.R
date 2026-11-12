# Autor: Dr. José Juan Torres Martínez
# Universidad Autónoma de San Luis Potosí (UASLP), México
#
# Figuras del artículo: Metagenomic characterization of Chlorella vulgaris and
# Scenedesmus acutus consortia (Algal Research, in preparation).
# Genera, a partir de los CSV de esta carpeta (valores tal como aparecen en
# la tesis de maestría de Abner Martínez-Constantino, FCQ-UASLP, agosto 2026):
#   Figura 1 — Comparación de los ensamblajes (contigs por umbral + N50)
#   Figura 2 — Bacterias más representativas por metagenoma (Kraken2)
#   Figura 3 — Tamaño de los bins que pasan el filtro (MAGs + nucleares)
#   Figura 4 — Evaluación BUSCO de los bins nucleares candidatos
#   Figura 5 — Anotación funcional (anotadas vs. hipotéticas + COG)
#   Figura 6 — Matrices de puntuación de funcionalidad (FS) por MAG
#   Figura 7 — Random Forest: importancia de variables + observados vs.
#              predichos (valores leídos de las figuras de la tesis)
#   Figura S0 — Identidad ANIm de BIN005 contra el panel de referencia
#
Sys.setlocale("LC_ALL", "C.UTF-8")

library(ggplot2)
library(dplyr)
library(readr)
library(scales)
library(tidyr)
library(patchwork)

# Carpeta del script: los CSV y las figuras viven junto a este archivo.
args <- commandArgs(trailingOnly = FALSE)
archivo <- grep("--file=", args, value = TRUE)
if (length(archivo) > 0) {
  carpeta <- dirname(sub("--file=", "", archivo[1]))
} else {
  carpeta <- getwd()
}
setwd(carpeta)

# Paleta Okabe-Ito (compatible con visión de color deficiente).
okabe <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999"
)
col_cv <- "#0072B2"
col_sa <- "#D55E00"

# Tema común: limpio, sin título dentro de la imagen, leyenda abajo,
# títulos de ejes en negrita.
tema_fig <- theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 11),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  )

# Etiquetas de taxón con género y especie en itálica (plotmath).
# "Pseudomonas oryzihabitans" -> italic('Pseudomonas oryzihabitans')
# "Microvirga sp. VF16"       -> italic('Microvirga')~sp.~'VF16'
# "Dyadobacter fermentans DSM 18053" -> italic('Dyadobacter fermentans')~'DSM 18053'
# La función devuelve una CADENA plotmath; etiquetas_taxon() la convierte
# en un vector de expresiones con un solo parse(text = ...) — un list de
# expresiones NO funciona como labels de una escala discreta en ggplot2.
etiqueta_taxon_str <- function(x) {
  partes <- strsplit(x, " ", fixed = TRUE)[[1]]
  genero <- partes[1]
  if (length(partes) == 1) {
    return(sprintf("italic('%s')", genero))
  }
  segundo <- partes[2]
  resto <- if (length(partes) > 2) paste(partes[3:length(partes)], collapse = " ") else ""
  if (segundo == "sp.") {
    expr <- sprintf("italic('%s')~sp.", genero)
    if (nzchar(resto)) expr <- paste0(expr, "~'", resto, "'")
    return(expr)
  }
  if (grepl("^[a-z]", segundo)) {
    expr <- sprintf("italic('%s %s')", genero, segundo)
    if (nzchar(resto)) expr <- paste0(expr, "~'", resto, "'")
    return(expr)
  }
  paste0(sprintf("italic('%s')", genero), "~'",
         paste(partes[2:length(partes)], collapse = " "), "'")
}
etiquetas_taxon <- function(v) parse(text = vapply(v, etiqueta_taxon_str, character(1)))

# Etiquetas de metagenoma con la especie en itálica (plotmath).
# Se usa en ejes discretos (Figuras 4 y 5A), en leyendas (Figuras 1, 5B
# y 7A-7B) y en los títulos de panel (Figuras 2, 3, 6 y 7C).
etiqueta_metagenoma <- function(v) {
  parse(text = vapply(v, function(x) {
    if (x == "C. vulgaris") "italic('C. vulgaris')"
    else if (x == "S. acutus") "italic('S. acutus')"
    else x
  }, character(1)))
}
labeller_metagenoma <- as_labeller(
  c("C. vulgaris" = "italic('C. vulgaris')", "S. acutus" = "italic('S. acutus')"),
  default = label_parsed)

# Etiquetas del panel ANI (Figura S0): cadena plotmath por referencia.
etiqueta_ani_str <- function(lbl) {
  if (grepl("reinhardtii", lbl)) return("italic('Chlamydomonas reinhardtii')~'CC-503 cw92'")
  if (grepl("tauri", lbl)) return("italic('Ostreococcus tauri')~'RCC4221'")
  if (grepl("LBCVB", lbl)) return("italic('C. vulgaris')~'LBCVB'~'(laboratory assembly)'")
  if (grepl("UTEX259", lbl)) return("italic('C. vulgaris')~'UTEX259'")
  if (grepl("NJ-7", lbl)) return("italic('C. vulgaris')~'NJ-7'")
  if (grepl("PKVL7422", lbl)) return("italic('C. vulgaris')~'PKVL7422'")
  if (grepl("GCA_211-11P", lbl)) return("italic('C. vulgaris')~'GCA_211-11P'")
  if (grepl("UTEX_395", lbl)) return("italic('C. vulgaris')~'UTEX_395'")
  lbl
}

# ---------------------------------------------------------------------------
# Figura 1 — Comparación de los ensamblajes
# ---------------------------------------------------------------------------
ens <- read_csv("datos_ensamblaje.csv", show_col_types = FALSE)
met <- read_csv("datos_metricas_ensamblaje.csv", show_col_types = FALSE)

ens$metagenome <- factor(ens$metagenome, levels = c("C. vulgaris", "S. acutus"))
ens$umbral <- factor(paste0(">=", ens$threshold_bp),
                     levels = c(">=0", ">=1000", ">=5000", ">=10000", ">=25000", ">=50000"))
p1a <- ggplot(ens, aes(x = umbral, y = contigs, color = metagenome, group = metagenome)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_y_continuous(labels = label_comma()) +
  scale_color_manual(values = c("C. vulgaris" = col_cv, "S. acutus" = col_sa),
                     labels = etiqueta_metagenoma) +
  labs(x = "Minimum contig length (bp)", y = "Number of contigs", color = "Metagenome") +
  tema_fig

n50 <- met %>% filter(metric %in% c("N50", "N90", "auN"))
n50$metric <- factor(n50$metric, levels = c("N50", "N90", "auN"))
n50$metagenome <- factor(n50$metagenome, levels = c("C. vulgaris", "S. acutus"))
p1b <- ggplot(n50, aes(x = metric, y = value, fill = metagenome)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_text(aes(label = comma(value)), position = position_dodge(width = 0.7),
            vjust = -0.4, size = 3.5) +
  scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.12))) +
  scale_fill_manual(values = c("C. vulgaris" = col_cv, "S. acutus" = col_sa),
                    labels = etiqueta_metagenoma) +
  labs(x = NULL, y = "Length (bp)", fill = "Metagenome") +
  tema_fig

fig1 <- p1a + p1b + plot_layout(ncol = 2, widths = c(1.25, 1))
ggsave("figura1_ensamblaje.png", fig1, width = 10, height = 4.5, dpi = 300, bg = "white")

# ---------------------------------------------------------------------------
# Figura 2 — Bacterias más representativas (Kraken2, lecturas)
# ---------------------------------------------------------------------------
tax <- read_csv("datos_taxonomia.csv", show_col_types = FALSE)
tax$metagenome <- factor(tax$metagenome, levels = c("C. vulgaris", "S. acutus"))
# reorder() es global: un taxón presente en ambos metagenomas
# (Ralstonia solanacearum) quedaría mal ordenado en cada panel.
# Se usa un factor con prefijo de metagenoma para ordenar dentro de cada panel.
tax <- tax %>%
  arrange(metagenome, reads) %>%
  mutate(taxon_f = factor(paste(metagenome, taxon, sep = "___"),
                          levels = paste(metagenome, taxon, sep = "___")))
etiquetas_taxon_prefijo <- function(v) {
  parse(text = vapply(v, function(lbl) {
    etiqueta_taxon_str(sub("^.*___", "", lbl))
  }, character(1)))
}

p2 <- ggplot(tax, aes(x = reads, y = taxon_f, fill = metagenome)) +
  geom_col() +
  facet_wrap(~metagenome, scales = "free_y", ncol = 1, labeller = labeller_metagenoma) +
  scale_x_continuous(trans = "log10", labels = label_comma()) +
  scale_y_discrete(labels = etiquetas_taxon_prefijo) +
  scale_fill_manual(values = c("C. vulgaris" = col_cv, "S. acutus" = col_sa)) +
  labs(x = "Reads (log scale)", y = NULL, fill = "Metagenome") +
  tema_fig +
  theme(legend.position = "none", axis.text.y = element_text(size = 10))
ggsave("figura2_taxonomia.png", p2, width = 9, height = 8, dpi = 300, bg = "white")

# ---------------------------------------------------------------------------
# Figura 3 — Tamaño de los bins que pasan el filtro
# ---------------------------------------------------------------------------
mags <- read_csv("datos_mags.csv", show_col_types = FALSE)
mags$metagenome <- factor(mags$metagenome, levels = c("C. vulgaris", "S. acutus"))
mags$tipo <- factor(ifelse(mags$type == "nuclear candidate",
                           "Nuclear candidate", "Bacterial"),
                    levels = c("Bacterial", "Nuclear candidate"))
mags <- mags %>%
  arrange(metagenome, size_Mb) %>%
  mutate(taxon_f = factor(paste(metagenome, taxon, sep = "___"),
                          levels = paste(metagenome, taxon, sep = "___")))

p3 <- ggplot(mags, aes(x = size_Mb, y = taxon_f, fill = tipo)) +
  geom_col() +
  geom_text(aes(label = paste0(bin, ": ", sprintf("%.3f", size_Mb))), hjust = -0.1, size = 3.2) +
  facet_wrap(~metagenome, scales = "free_y", ncol = 1, labeller = labeller_metagenoma) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.28))) +
  scale_y_discrete(labels = etiquetas_taxon_prefijo) +
  scale_fill_manual(values = c("Bacterial" = "#56B4E9", "Nuclear candidate" = "#D55E00")) +
  labs(x = "Bin size (Mb)", y = NULL, fill = "Bin type") +
  tema_fig +
  theme(axis.text.y = element_text(size = 10))
ggsave("figura3_mags.png", p3, width = 10, height = 8, dpi = 300, bg = "white")

# ---------------------------------------------------------------------------
# Figura 4 — BUSCO de los bins nucleares candidatos
# ---------------------------------------------------------------------------
busco <- read_csv("datos_busco.csv", show_col_types = FALSE)
busco$metagenome <- factor(busco$metagenome, levels = c("C. vulgaris", "S. acutus"))
busco$category <- factor(busco$category,
  levels = c("Complete (single-copy)", "Duplicated", "Fragmented", "Missing"))
p4 <- ggplot(busco, aes(x = metagenome, y = percent, fill = category)) +
  geom_col() +
  geom_text(data = function(d) d[d$percent >= 5, ],
            aes(label = sprintf("%.1f%%", percent)),
            position = position_stack(vjust = 0.5), size = 3.5, color = "white",
            fontface = "bold") +
  scale_fill_manual(values = c("Complete (single-copy)" = "#009E73",
                                "Duplicated" = "#56B4E9",
                                "Fragmented" = "#E69F00",
                                "Missing" = "#999999")) +
  scale_x_discrete(labels = etiqueta_metagenoma) +
  labs(x = NULL, y = "BUSCO genes (%)", fill = "BUSCO category") +
  tema_fig
ggsave("figura4_busco.png", p4, width = 7, height = 5, dpi = 300, bg = "white")

# ---------------------------------------------------------------------------
# Figura 5 — Anotación funcional (A) y COG (B)
# ---------------------------------------------------------------------------
anot <- read_csv("datos_anotacion.csv", show_col_types = FALSE)
anot$metagenome <- factor(anot$metagenome, levels = c("C. vulgaris", "S. acutus"))
anot$category <- factor(anot$category,
  levels = c("Annotated proteins", "Hypothetical proteins"))
p5a <- ggplot(anot, aes(x = metagenome, y = count, fill = category)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_text(aes(label = comma(count)), position = position_dodge(width = 0.7),
            vjust = -0.4, size = 3.2) +
  scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.15))) +
  scale_x_discrete(labels = etiqueta_metagenoma) +
  scale_fill_manual(values = c("Annotated proteins" = "#009E73",
                                "Hypothetical proteins" = "#999999")) +
  labs(x = NULL, y = "Proteins", fill = NULL) +
  tema_fig

cog <- read_csv("datos_cog.csv", show_col_types = FALSE)
cog <- cog %>% filter(reported == "yes")
cog$metagenome <- factor(cog$metagenome, levels = c("C. vulgaris", "S. acutus"))
p5b <- ggplot(cog, aes(x = cog_category, y = percent, fill = metagenome)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_text(aes(label = sprintf("%.2f%%", percent)),
            position = position_dodge(width = 0.7), vjust = -0.4, size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = c("C. vulgaris" = col_cv, "S. acutus" = col_sa),
                    labels = etiqueta_metagenoma) +
  labs(x = NULL, y = "Proteins (%)", fill = "Metagenome") +
  tema_fig +
  theme(axis.text.x = element_text(size = 10))

fig5 <- p5a + p5b + plot_layout(ncol = 2, widths = c(1, 1))
ggsave("figura5_anotacion_cog.png", fig5, width = 11, height = 4.5, dpi = 300, bg = "white")

# ---------------------------------------------------------------------------
# Figura 6 — Matrices de puntuación de funcionalidad (FS)
# ---------------------------------------------------------------------------
fs <- read_csv("datos_fs.csv", show_col_types = FALSE)
fs$metagenome <- factor(fs$metagenome, levels = c("C. vulgaris", "S. acutus"))
fs_largo <- fs %>%
  pivot_longer(cols = c(biotin, riboflavin, quorum_sensing, siderophores),
               names_to = "pathway", values_to = "fs")
fs_largo$pathway <- factor(fs_largo$pathway,
  levels = c("biotin", "riboflavin", "quorum_sensing", "siderophores"),
  labels = c("Biotin", "Riboflavin", "Quorum sensing", "Siderophores"))

p6 <- ggplot(fs_largo, aes(x = pathway, y = mag, fill = fs)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", fs)), size = 3) +
  facet_wrap(~metagenome, ncol = 1, scales = "free_y", labeller = labeller_metagenoma) +
  scale_fill_gradient(low = "#F0F0F0", high = "#0072B2", limits = c(0, 1),
                      name = "Functional score") +
  labs(x = NULL, y = "MAG") +
  tema_fig +
  theme(axis.text.x = element_text(size = 11, face = "bold"))
ggsave("figura6_fs.png", p6, width = 8, height = 9, dpi = 300, bg = "white")

# ---------------------------------------------------------------------------
# Figura 7 — Random Forest (valores APROXIMADOS leídos de las figuras)
# ---------------------------------------------------------------------------
rf <- read_csv("datos_rf_importancia.csv", show_col_types = FALSE)
rf$metagenome <- factor(rf$metagenome, levels = c("C. vulgaris", "S. acutus"))
rf$pathway <- factor(rf$pathway,
  levels = c("Biotin", "Riboflavin", "Quorum sensing", "Siderophores"))

p7a <- ggplot(rf, aes(x = pathway, y = inc_mse, fill = metagenome)) +
  scale_x_discrete(labels = c("Biotin", "Riboflavin", "Quorum\nsensing", "Siderophores")) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  geom_text(aes(label = inc_mse), position = position_dodge(width = 0.7),
            vjust = -0.4, size = 3.2) +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  scale_fill_manual(values = c("C. vulgaris" = col_cv, "S. acutus" = col_sa),
                    labels = etiqueta_metagenoma) +
  labs(x = NULL, y = "%IncMSE", fill = "Metagenome") +
  tema_fig

p7b <- ggplot(rf, aes(x = pathway, y = inc_node_purity, fill = metagenome)) +
  scale_x_discrete(labels = c("Biotin", "Riboflavin", "Quorum\nsensing", "Siderophores")) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_text(aes(label = sprintf("%.2f", inc_node_purity)),
            position = position_dodge(width = 0.7), vjust = -0.4, size = 3.2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = c("C. vulgaris" = col_cv, "S. acutus" = col_sa),
                    labels = etiqueta_metagenoma) +
  labs(x = NULL, y = "IncNodePurity", fill = "Metagenome") +
  tema_fig

pred <- read_csv("datos_rf_predicho.csv", show_col_types = FALSE)
pred$metagenome <- factor(pred$metagenome, levels = c("C. vulgaris", "S. acutus"))
p7c <- ggplot(pred, aes(x = (observed_low + observed_high) / 2,
                        y = (predicted_low + predicted_high) / 2)) +
  geom_errorbar(aes(ymin = predicted_low, ymax = predicted_high),
                width = 0.03, color = "#0072B2") +
  geom_errorbarh(aes(xmin = observed_low, xmax = observed_high),
                  height = 0.03, color = "#0072B2") +
  geom_point(size = 3, color = "#0072B2") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#999999") +
  facet_wrap(~metagenome, ncol = 2, labeller = labeller_metagenoma) +
  coord_cartesian(xlim = c(0, 1.05), ylim = c(0, 1.05)) +
  labs(x = "Observed influence score", y = "Predicted influence score") +
  tema_fig

fig7 <- (p7a + p7b) / p7c + plot_layout(heights = c(1, 1))
ggsave("figura7_random_forest.png", fig7, width = 10, height = 8, dpi = 300, bg = "white")

# ---------------------------------------------------------------------------
# Figura S0 — Identidad ANIm de BIN005 (leída de la Fig. 12 de la tesis)
# ---------------------------------------------------------------------------
ani <- read_csv("datos_ani.csv", show_col_types = FALSE)
ani <- ani %>% arrange(desc(ani_percent))
pS0 <- ggplot(ani, aes(x = ani_percent,
                       y = reorder(reference, ani_percent),
                       fill = ani_percent >= 95)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.2f%%", ani_percent)), hjust = -0.1, size = 3.2) +
  geom_vline(xintercept = 95, linetype = "dashed", color = "#D55E00") +
  scale_x_continuous(limits = c(0, 108), expand = expansion(mult = c(0, 0.02))) +
  scale_y_discrete(labels = function(v) parse(text = vapply(v, etiqueta_ani_str, character(1)))) +
  scale_fill_manual(values = c("TRUE" = "#009E73", "FALSE" = "#999999")) +
  labs(x = "ANIm identity (%)", y = NULL, fill = ">= 95%") +
  tema_fig +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 10))
ggsave("figuraS0_ani.png", pS0, width = 9, height = 5.5, dpi = 300, bg = "white")

message("Figuras generadas: figura1_ensamblaje.png, figura2_taxonomia.png, ",
        "figura3_mags.png, figura4_busco.png, figura5_anotacion_cog.png, ",
        "figura6_fs.png, figura7_random_forest.png, figuraS0_ani.png")
