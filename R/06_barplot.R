# Make barplots of
# genes : L10RA, CD160, FABP6, TLR6, CD274 (PD-L1), KYNU, ASCL2, SERPINI1, MSX2, sirt2, C2, C1R, BEST4
# Samples: Postbiotics and PBS
library("tidyr")
library("dplyr")
library("ggplot2")
library("ggpubr")
library("rstatix")
res <- readRDS("output/limma_juanjo.RDS")
counts <- readRDS("output/counts.RDS")
comp <- readRDS("output/samples_comparisons.RDS")
meta <- readRDS("output/pheno.RDS")

# * Create matrix for significance ####
diff <- matrix(dimnames = dimnames(res$fc), ncol = ncol(res$fc), nrow = nrow(res$fc))
diff[res$p < 0.05 & res$fc > 0] <- "UP"
diff[res$p < 0.05 & res$fc < 0] <- "DW"
diff[res$fdr < 0.05 & res$fc > 0] <- "UUP"
diff[res$fdr < 0.05 & res$fc < 0] <- "DDW"
colnames(diff) <- gsub("fc_", "", colnames(diff))

# Select samples
samples <- meta$`Macrogen SAMPLE NAME`[(!is.na(meta$PB) & is.na(meta$estimul)) | meta$cond == "PBS"]

# Normalize
cn <- integration::norm_RNAseq(counts[, samples])

# Select genes
genes <- rownames(cn)
g <- strsplit(genes, "\\.[0-9]*_")
g <- t(simplify2array(g))
colnames(g) <- c("ENSEMBL", "name")
g <- cbind(genes = genes, g)
names_genes <- c("IL10RA", "CD160", "FABP6", "TLR6", "CD274", "KYNU",
                 "ASCL2", "SERPINI1", "MSX2", "SIRT2", "C2", "C1R", "BEST4", "ORC1", "RIMBP3")
sel <- g[g[, "name"] %in% names_genes, ]

counts_df <- as.data.frame(cn[ sel[, 1], samples]) %>%
  tibble::rownames_to_column() %>%
  pivot_longer(cols = -rowname) %>%
  left_join(meta, by = c("name" = "Macrogen SAMPLE NAME")) %>%
  left_join(as.data.frame(sel), by = c("rowname" = "genes"))

ggplot(counts_df) +
  geom_abline(slope = 0, intercept = 0, col = "black") +
  geom_jitter(aes(cond, value, col = name.y, group = cond), height = 0) +
  facet_wrap(~name.y) +
  labs(col = "Gene", y = "Normalized expression", x = "Condition") +
  theme_minimal()
ggsave("Figures/points_genes.png")
my_comparisons <- list( c(1, 2), c(1, 3), c(1, 4))
a <- counts_df %>%
  # group_by(cond, name.y) %>%
  # summarize(mean_se(value)) %>%
  # ungroup() %>%
  mutate(cond = if_else(cond == "PBS", "vehicle", cond)) %>%
  mutate(cond = forcats::fct_relevel(cond, c("vehicle", "PBK12", "PBNissle", "PBSth"))) %>%
  # filter(name.y == "CD160") %>%
  ggboxplot(x = "cond", y = "value", col = "name.y", add = "jitter") +
  # geom_point(aes(col = name.y)) +
  # geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.5) +
  facet_wrap(~name.y, scales = "free") +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.adj") # Add pairwise comparisons p-value
  labs(fill = "Gene", y = "Normalized expression", x = element_blank()) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave("Figures/bars_genes_signif.png")


col <- grep("PBS_vs|_PBS$", colnames(res$p), value = TRUE)[c(10, 7:9)]
p_df <- res$p[sel[, 1], col] %>%
  as_tibble(rownames = "rowname") %>%
  left_join(as.data.frame(sel), by = c("rowname" = "genes"))
real_p <- p_df %>%
  select(-rowname, -ENSEMBL) %>%
  pivot_longer(cols = -name, names_to = "comparative", values_to = "p") %>%
  mutate(group1 = gsub("p_(.+)_vs_.+", "\\1", comparative),
         group2 = gsub("p_(.+)_vs_(.+)", "\\2", comparative)) %>%
  filter(group2 != "Control")


stat.test <- counts_df %>%
  group_by(name.y) %>%
  t_test(value~cond, ref.group = "PBS") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  add_y_position() %>%
  full_join(real_p, by  = c("group1" = "group2", "group2" = "group1", "name.y" = "name"))

stat.test2 <- stat.test %>%
  select(-p.x, -p.adj) %>%
  rename(p = p.y) %>%
  group_by(name.y) %>%
  mutate(p.adj = p.adjust(p, method = "bonferroni")) %>%
  ungroup() %>%
  mutate(group1 = if_else(group1 == "PBS", "vehicle", group1)) %>%
  mutate(p.formatted = scales::pvalue(p, accuracy= 0.001, add_p = TRUE),
         p.adj.formatted = scales::pvalue(p.adj, accuracy= 0.001, add_p = TRUE))

counts_df %>%
  mutate(cond = if_else(cond == "PBS", "vehicle", cond)) %>%
  mutate(cond = forcats::fct_relevel(cond, c("vehicle", "PBK12", "PBNissle", "PBSth"))) %>%
  ggboxplot(x = "cond", y = "value", facet.by = "name.y", color = "name.y", scales = "free") +
  rotate_x_text(angle = 60) +
  stat_pvalue_manual(stat.test2, label = "p.formatted", vjust = 1.5) +
  labs(fill = "Gene", y = "Normalized expression", x = element_blank()) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave("Figures/bars_genes_signif_macros.png")

counts_df %>%
  mutate(cond = if_else(cond == "PBS", "vehicle", cond)) %>%
  mutate(cond = forcats::fct_relevel(cond, c("vehicle", "PBK12", "PBNissle", "PBSth"))) %>%
  select(name.y, value, rowname, cond, SAMPLE) %>%
  writexl::write_xlsx("output/valors_bars_genes_signif_macros.xlsx")
