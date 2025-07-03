############DDR#######################

# ----------------- Load Libraries -----------------
library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)

# ----------------- DDR Entrez IDs -----------------
entrez_ids <- c(
  10111, 1017, 1019, 1020, 1021, 1026, 1027, 10912, 11011, 1111,
  11200, 1385, 1643, 1647, 1869, 207, 2177, 25, 27113, 27244,
  3014, 317, 355, 4193, 4292, 4361, 4609, 4616, 4683, 472, 50484,
  5366, 5371, 54205, 545, 55367, 5591, 581, 5810, 5883, 5884,
  5888, 5893, 5925, 595, 5981, 6118, 637, 672, 7157, 7799,
  8243, 836, 841, 84126, 842, 8795, 891, 894, 896, 898,
  9133, 9134, 983, 9874, 993, 995, 5916
)

# ----------------- Load CorMotif Groups -----------------
# 0.1 µM
prob_1_0.1 <- as.character(read.csv("data/prob_1_0.1.csv")$Entrez_ID)
prob_2_0.1 <- as.character(read.csv("data/prob_2_0.1.csv")$Entrez_ID)
prob_3_0.1 <- as.character(read.csv("data/prob_3_0.1.csv")$Entrez_ID)

# 0.5 µM
prob_1_0.5 <- as.character(read.csv("data/prob_1_0.5.csv")$Entrez_ID)
prob_2_0.5 <- as.character(read.csv("data/prob_2_0.5.csv")$Entrez_ID)
prob_3_0.5 <- as.character(read.csv("data/prob_3_0.5.csv")$Entrez_ID)
prob_4_0.5 <- as.character(read.csv("data/prob_4_0.5.csv")$Entrez_ID)
prob_5_0.5 <- as.character(read.csv("data/prob_5_0.5.csv")$Entrez_ID)

# ----------------- Annotate Response Groups -----------------
df_0.1 <- data.frame(Entrez_ID = unique(c(prob_1_0.1, prob_2_0.1, prob_3_0.1))) %>%
  filter(Entrez_ID %in% as.character(entrez_ids)) %>%
  mutate(
    Response_Group = case_when(
      Entrez_ID %in% prob_1_0.1 ~ "Non response",
      Entrez_ID %in% prob_2_0.1 ~ "CX-DOX mid-late response",
      Entrez_ID %in% prob_3_0.1 ~ "DOX only mid-late response"
    ),
    Concentration = "0.1 µM"
  )

df_0.5 <- data.frame(Entrez_ID = unique(c(prob_1_0.5, prob_2_0.5, prob_3_0.5, prob_4_0.5, prob_5_0.5))) %>%
  filter(Entrez_ID %in% as.character(entrez_ids)) %>%
  mutate(
    Response_Group = case_when(
      Entrez_ID %in% prob_1_0.5 ~ "Non response",
      Entrez_ID %in% prob_2_0.5 ~ "DOX-specific response",
      Entrez_ID %in% prob_3_0.5 ~ "DOX only mid-late response",
      Entrez_ID %in% prob_4_0.5 ~ "CX total + DOX early response",
      Entrez_ID %in% prob_5_0.5 ~ "DOX early + CX-DOX mid-late response"
    ),
    Concentration = "0.5 µM"
  )

df_combined <- bind_rows(df_0.1, df_0.5)

# ----------------- Load and Aggregate DEG Data -----------------
deg_files <- list(
  "CX_0.1_3" = "data/DEGs/Toptable_CX_0.1_3.csv",
  "CX_0.1_24" = "data/DEGs/Toptable_CX_0.1_24.csv",
  "CX_0.1_48" = "data/DEGs/Toptable_CX_0.1_48.csv",
  "CX_0.5_3" = "data/DEGs/Toptable_CX_0.5_3.csv",
  "CX_0.5_24" = "data/DEGs/Toptable_CX_0.5_24.csv",
  "CX_0.5_48" = "data/DEGs/Toptable_CX_0.5_48.csv",
  "DOX_0.1_3" = "data/DEGs/Toptable_DOX_0.1_3.csv",
  "DOX_0.1_24" = "data/DEGs/Toptable_DOX_0.1_24.csv",
  "DOX_0.1_48" = "data/DEGs/Toptable_DOX_0.1_48.csv",
  "DOX_0.5_3" = "data/DEGs/Toptable_DOX_0.5_3.csv",
  "DOX_0.5_24" = "data/DEGs/Toptable_DOX_0.5_24.csv",
  "DOX_0.5_48" = "data/DEGs/Toptable_DOX_0.5_48.csv"
)

deg_all <- map2_dfr(deg_files, names(deg_files), ~{
  read_csv(.x, show_col_types = FALSE) %>%
    mutate(Condition = .y)
}) %>%
  mutate(
    Entrez_ID = as.character(Entrez_ID),
    Concentration = case_when(
      str_detect(Condition, "0.1") ~ "0.1 µM",
      str_detect(Condition, "0.5") ~ "0.5 µM",
      TRUE ~ NA_character_
    )
  )

deg_mean <- deg_all %>%
  group_by(Entrez_ID, Concentration) %>%
  summarise(
    logFC = mean(logFC, na.rm = TRUE),
    .groups = "drop"
  )

# ----------------- Merge with CorMotif DDR Genes -----------------
ddr_logfc <- inner_join(df_combined, deg_mean, by = c("Entrez_ID", "Concentration")) %>%
  mutate(
    Gene = mapIds(org.Hs.eg.db, keys = Entrez_ID,
                  column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"),
    Group = paste(Response_Group, Concentration)
  )

# ----------------- Define order -----------------
group_levels <- c(
  "Non response",
  "CX-DOX mid-late response",
  "DOX-specific response",
  "DOX only mid-late response",
  "CX total + DOX early response",
  "DOX early + CX-DOX mid-late response"
)

ddr_logfc$Response_Group <- factor(ddr_logfc$Response_Group, levels = group_levels)

# ----------------- Wilcoxon Test -----------------
comparison_map <- list(
  "CX-DOX mid-late response 0.1 µM" = "Non response 0.1 µM",
  "DOX only mid-late response 0.1 µM" = "Non response 0.1 µM",
  "DOX-specific response 0.5 µM" = "Non response 0.5 µM",
  "DOX only mid-late response 0.5 µM" = "Non response 0.5 µM",
  "CX total + DOX early response 0.5 µM" = "Non response 0.5 µM",
  "DOX early + CX-DOX mid-late response 0.5 µM" = "Non response 0.5 µM"
)

star_df <- lapply(names(comparison_map), function(resp_group) {
  control_group <- comparison_map[[resp_group]]
  resp_vals <- ddr_logfc$logFC[ddr_logfc$Group == resp_group]
  control_vals <- ddr_logfc$logFC[ddr_logfc$Group == control_group]

  if (length(resp_vals) >= 2 && length(control_vals) >= 2) {
    test_result <- wilcox.test(resp_vals, control_vals)
    pval <- test_result$p.value
    if (pval < 0.05) {
      label <- case_when(
        pval < 0.001 ~ "***",
        pval < 0.01 ~ "**",
        TRUE ~ "*"
      )
      y_pos <- max(resp_vals, na.rm = TRUE) + 0.5
      return(data.frame(Group = resp_group, y_pos = y_pos, label = label, P_Value = signif(pval, 4)))
    }
  }
  return(NULL)
}) %>% bind_rows()

# ----------------- Prepare Annotation -----------------
label_data <- star_df %>%
  separate(Group, into = c("Response_Group", "Concentration"), sep = " (?=0\\.)") %>%
  mutate(Response_Group = factor(Response_Group, levels = group_levels))

# ----------------- Plot -----------------
ggplot(ddr_logfc, aes(x = Response_Group, y = logFC, fill = Response_Group)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7) +
  geom_text(data = label_data, aes(x = Response_Group, y = y_pos, label = label),
            inherit.aes = FALSE, size = 5, color = "black") +
  facet_wrap(~Concentration, scales = "free_x") +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Mean logFC of DDR Genes\nAcross CorMotif Response Groups",
    x = "Response Group",
    y = "Mean Log Fold Change",
    fill = "Response Group"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 13),
    strip.text = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)
  )







############ P53 TARGET GENE ANALYSIS #######################

# ----------------- Load Libraries -----------------
library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)

# ----------------- P53 Entrez IDs -----------------
entrez_ids <- c(
  1026, 50484, 4193, 9766, 9518, 7832, 1643, 1647, 1263, 57103, 51065, 8795, 51499, 64393, 581, 5228, 5429, 8493, 55959,
  7508, 64782, 282991, 355, 53836, 4814, 10769, 9050, 27244, 9540, 94241, 26154, 57763, 900, 26999, 55332, 26263, 23479,
  23612, 29950, 9618, 10346, 8824, 134147, 55294, 22824, 4254, 6560, 467, 27113, 60492, 8444, 60401, 1969, 220965, 2232,
  3976, 55191, 84284, 93129, 5564, 7803, 83667, 7779, 132671, 7039, 51768, 137695, 93134, 7633, 10973, 340485, 307,
  27350, 23245, 3732, 29965, 1363, 1435, 196513, 8507, 8061, 2517, 51278, 53354, 54858, 23228, 5366, 5912, 6236, 51222,
  26152, 59, 1907, 50650, 91012, 780, 9249, 11072, 144455, 64787, 116151, 27165, 2876, 57822, 55733, 57722, 121457,
  375449, 85377, 4851, 5875, 127544, 29901, 84958, 8797, 8793, 441631, 220001, 54541, 5889, 5054, 25816, 25987, 5111,
  98, 317, 598, 604, 10904, 1294, 80315, 53944, 1606, 2770, 3628, 3675, 3985, 4035, 4163, 84552, 29085, 55367, 5371,
  5791, 54884, 5980, 8794, 1462, 50808, 220, 583, 694, 1056, 9076, 10978, 54677, 1612, 55040, 114907, 2274, 127707,
  4000, 8079, 4646, 4747, 27445, 5143, 80055, 79156, 5360, 5364, 23654, 5565, 5613, 5625, 10076, 56963, 6004, 390,
  255488, 6326, 6330, 23513, 7869, 283130, 204962, 83959, 6548, 6774, 9263, 10228, 22954, 10475, 85363, 494514, 10142,
  79714, 1006, 8446, 9648, 79828, 5507, 55240, 63874, 25841, 9289, 84883, 154810, 51321, 421, 8553, 655, 119032, 84280,
  10950, 824, 839, 57828, 857, 8812, 8837, 94027, 113189, 22837, 132864, 10898, 3300, 81704, 1847, 1849, 1947, 9538,
  24139, 5168, 147965, 115548, 9873, 23768, 2632, 2817, 3280, 3265, 23308, 3490, 51477, 182, 3856, 8844, 144811, 9404,
  4043, 9848, 2872, 23041, 740, 343263, 4638, 26509, 4792, 22861, 57523, 55214, 80025, 164091, 57060, 64065, 51090,
  5453, 8496, 333926, 55671, 5900, 55544, 23179, 8601, 389, 6223, 55800, 6385, 4088, 6643, 122809, 257397, 285343,
  7011, 54790, 374618, 55362, 51754, 7157, 9537, 22906, 7205, 80705, 219699, 55245, 83719, 7748, 25946, 118738
)

# ----------------- Load CorMotif Groups -----------------
# 0.1 µM
prob_1_0.1 <- as.character(read.csv("data/prob_1_0.1.csv")$Entrez_ID)
prob_2_0.1 <- as.character(read.csv("data/prob_2_0.1.csv")$Entrez_ID)
prob_3_0.1 <- as.character(read.csv("data/prob_3_0.1.csv")$Entrez_ID)

# 0.5 µM
prob_1_0.5 <- as.character(read.csv("data/prob_1_0.5.csv")$Entrez_ID)
prob_2_0.5 <- as.character(read.csv("data/prob_2_0.5.csv")$Entrez_ID)
prob_3_0.5 <- as.character(read.csv("data/prob_3_0.5.csv")$Entrez_ID)
prob_4_0.5 <- as.character(read.csv("data/prob_4_0.5.csv")$Entrez_ID)
prob_5_0.5 <- as.character(read.csv("data/prob_5_0.5.csv")$Entrez_ID)

# ----------------- Annotate Response Groups -----------------
df_0.1 <- data.frame(Entrez_ID = unique(c(prob_1_0.1, prob_2_0.1, prob_3_0.1))) %>%
  filter(Entrez_ID %in% as.character(entrez_ids)) %>%
  mutate(
    Response_Group = case_when(
      Entrez_ID %in% prob_1_0.1 ~ "Non response",
      Entrez_ID %in% prob_2_0.1 ~ "CX-DOX mid-late response",
      Entrez_ID %in% prob_3_0.1 ~ "DOX only mid-late response"
    ),
    Concentration = "0.1 µM"
  )

df_0.5 <- data.frame(Entrez_ID = unique(c(prob_1_0.5, prob_2_0.5, prob_3_0.5, prob_4_0.5, prob_5_0.5))) %>%
  filter(Entrez_ID %in% as.character(entrez_ids)) %>%
  mutate(
    Response_Group = case_when(
      Entrez_ID %in% prob_1_0.5 ~ "Non response",
      Entrez_ID %in% prob_2_0.5 ~ "DOX-specific response",
      Entrez_ID %in% prob_3_0.5 ~ "DOX only mid-late response",
      Entrez_ID %in% prob_4_0.5 ~ "CX total + DOX early response",
      Entrez_ID %in% prob_5_0.5 ~ "DOX early + CX-DOX mid-late response"
    ),
    Concentration = "0.5 µM"
  )

df_combined <- bind_rows(df_0.1, df_0.5)

# ----------------- Load and Aggregate DEG Data -----------------
deg_files <- list(
  "CX_0.1_3" = "data/DEGs/Toptable_CX_0.1_3.csv",
  "CX_0.1_24" = "data/DEGs/Toptable_CX_0.1_24.csv",
  "CX_0.1_48" = "data/DEGs/Toptable_CX_0.1_48.csv",
  "CX_0.5_3" = "data/DEGs/Toptable_CX_0.5_3.csv",
  "CX_0.5_24" = "data/DEGs/Toptable_CX_0.5_24.csv",
  "CX_0.5_48" = "data/DEGs/Toptable_CX_0.5_48.csv",
  "DOX_0.1_3" = "data/DEGs/Toptable_DOX_0.1_3.csv",
  "DOX_0.1_24" = "data/DEGs/Toptable_DOX_0.1_24.csv",
  "DOX_0.1_48" = "data/DEGs/Toptable_DOX_0.1_48.csv",
  "DOX_0.5_3" = "data/DEGs/Toptable_DOX_0.5_3.csv",
  "DOX_0.5_24" = "data/DEGs/Toptable_DOX_0.5_24.csv",
  "DOX_0.5_48" = "data/DEGs/Toptable_DOX_0.5_48.csv"
)

deg_all <- map2_dfr(deg_files, names(deg_files), ~{
  read_csv(.x, show_col_types = FALSE) %>%
    mutate(Condition = .y)
}) %>%
  mutate(
    Entrez_ID = as.character(Entrez_ID),
    Concentration = case_when(
      str_detect(Condition, "0.1") ~ "0.1 µM",
      str_detect(Condition, "0.5") ~ "0.5 µM",
      TRUE ~ NA_character_
    )
  )

deg_mean <- deg_all %>%
  group_by(Entrez_ID, Concentration) %>%
  summarise(logFC = mean(logFC, na.rm = TRUE), .groups = "drop")

# ----------------- Merge with CorMotif P53 Genes -----------------
p53_logfc <- inner_join(df_combined, deg_mean, by = c("Entrez_ID", "Concentration")) %>%
  mutate(
    Gene = mapIds(org.Hs.eg.db, keys = Entrez_ID,
                  column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"),
    Group = paste(Response_Group, Concentration)
  )

# ----------------- Define order -----------------
group_levels <- c(
  "Non response",
  "CX-DOX mid-late response",
  "DOX-specific response",
  "DOX only mid-late response",
  "CX total + DOX early response",
  "DOX early + CX-DOX mid-late response"
)

p53_logfc$Response_Group <- factor(p53_logfc$Response_Group, levels = group_levels)

# ----------------- Wilcoxon Test -----------------
comparison_map <- list(
  "CX-DOX mid-late response 0.1 µM" = "Non response 0.1 µM",
  "DOX only mid-late response 0.1 µM" = "Non response 0.1 µM",
  "DOX-specific response 0.5 µM" = "Non response 0.5 µM",
  "DOX only mid-late response 0.5 µM" = "Non response 0.5 µM",
  "CX total + DOX early response 0.5 µM" = "Non response 0.5 µM",
  "DOX early + CX-DOX mid-late response 0.5 µM" = "Non response 0.5 µM"
)

star_df <- lapply(names(comparison_map), function(resp_group) {
  control_group <- comparison_map[[resp_group]]
  resp_vals <- p53_logfc$logFC[p53_logfc$Group == resp_group]
  control_vals <- p53_logfc$logFC[p53_logfc$Group == control_group]

  if (length(resp_vals) >= 2 && length(control_vals) >= 2) {
    test_result <- wilcox.test(resp_vals, control_vals)
    pval <- test_result$p.value
    if (pval < 0.05) {
      label <- case_when(
        pval < 0.001 ~ "***",
        pval < 0.01 ~ "**",
        TRUE ~ "*"
      )
      y_pos <- max(resp_vals, na.rm = TRUE) + 0.5
      return(data.frame(Group = resp_group, y_pos = y_pos, label = label, P_Value = signif(pval, 4)))
    }
  }
  return(NULL)
}) %>% bind_rows()


# ----------------- Prepare Annotation -----------------
label_data <- star_df %>%
  separate(Group, into = c("Response_Group", "Concentration"), sep = " (?=0\\.)") %>%
  mutate(Response_Group = factor(Response_Group, levels = group_levels))




# ----------------- Plot -----------------
ggplot(p53_logfc, aes(x = Response_Group, y = logFC, fill = Response_Group)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7) +
  geom_text(data = label_data, aes(x = Response_Group, y = y_pos, label = label),
            inherit.aes = FALSE, size = 5, color = "black") +
  facet_wrap(~Concentration, scales = "free_x") +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Mean logFC of P53 Target Genes\nAcross CorMotif Response Groups",
    x = "Response Group",
    y = "Mean Log Fold Change",
    fill = "Response Group"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 13),
    strip.text = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)
  )







############ TOP2B TARGET GENE ANALYSIS #######################

# ----------------- Load Libraries -----------------
library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)

# ----------------- TOP2B Entrez IDs -----------------
top2b_entrez_ids <- read_csv("data/TOP2B_target_mapped.csv", show_col_types = FALSE) %>%
  pull(Entrez_ID) %>%
  unique() %>%
  as.character()

# ----------------- Load CorMotif Groups -----------------
# 0.1 µM
prob_1_0.1 <- as.character(read.csv("data/prob_1_0.1.csv")$Entrez_ID)
prob_2_0.1 <- as.character(read.csv("data/prob_2_0.1.csv")$Entrez_ID)
prob_3_0.1 <- as.character(read.csv("data/prob_3_0.1.csv")$Entrez_ID)

# 0.5 µM
prob_1_0.5 <- as.character(read.csv("data/prob_1_0.5.csv")$Entrez_ID)
prob_2_0.5 <- as.character(read.csv("data/prob_2_0.5.csv")$Entrez_ID)
prob_3_0.5 <- as.character(read.csv("data/prob_3_0.5.csv")$Entrez_ID)
prob_4_0.5 <- as.character(read.csv("data/prob_4_0.5.csv")$Entrez_ID)
prob_5_0.5 <- as.character(read.csv("data/prob_5_0.5.csv")$Entrez_ID)

# ----------------- Annotate Response Groups -----------------
df_0.1 <- data.frame(Entrez_ID = unique(c(prob_1_0.1, prob_2_0.1, prob_3_0.1))) %>%
  filter(Entrez_ID %in% top2b_entrez_ids) %>%
  mutate(
    Response_Group = case_when(
      Entrez_ID %in% prob_1_0.1 ~ "Non response",
      Entrez_ID %in% prob_2_0.1 ~ "CX-DOX mid-late response",
      Entrez_ID %in% prob_3_0.1 ~ "DOX only mid-late response"
    ),
    Concentration = "0.1 µM"
  )

df_0.5 <- data.frame(Entrez_ID = unique(c(prob_1_0.5, prob_2_0.5, prob_3_0.5, prob_4_0.5, prob_5_0.5))) %>%
  filter(Entrez_ID %in% top2b_entrez_ids) %>%
  mutate(
    Response_Group = case_when(
      Entrez_ID %in% prob_1_0.5 ~ "Non response",
      Entrez_ID %in% prob_2_0.5 ~ "DOX-specific response",
      Entrez_ID %in% prob_3_0.5 ~ "DOX only mid-late response",
      Entrez_ID %in% prob_4_0.5 ~ "CX total + DOX early response",
      Entrez_ID %in% prob_5_0.5 ~ "DOX early + CX-DOX mid-late response"
    ),
    Concentration = "0.5 µM"
  )

df_combined <- bind_rows(df_0.1, df_0.5)

# ----------------- Load and Aggregate DEG Data -----------------
deg_files <- list(
  "CX_0.1_3" = "data/DEGs/Toptable_CX_0.1_3.csv",
  "CX_0.1_24" = "data/DEGs/Toptable_CX_0.1_24.csv",
  "CX_0.1_48" = "data/DEGs/Toptable_CX_0.1_48.csv",
  "CX_0.5_3" = "data/DEGs/Toptable_CX_0.5_3.csv",
  "CX_0.5_24" = "data/DEGs/Toptable_CX_0.5_24.csv",
  "CX_0.5_48" = "data/DEGs/Toptable_CX_0.5_48.csv",
  "DOX_0.1_3" = "data/DEGs/Toptable_DOX_0.1_3.csv",
  "DOX_0.1_24" = "data/DEGs/Toptable_DOX_0.1_24.csv",
  "DOX_0.1_48" = "data/DEGs/Toptable_DOX_0.1_48.csv",
  "DOX_0.5_3" = "data/DEGs/Toptable_DOX_0.5_3.csv",
  "DOX_0.5_24" = "data/DEGs/Toptable_DOX_0.5_24.csv",
  "DOX_0.5_48" = "data/DEGs/Toptable_DOX_0.5_48.csv"
)

deg_all <- map2_dfr(deg_files, names(deg_files), ~{
  read_csv(.x, show_col_types = FALSE) %>%
    mutate(Condition = .y)
}) %>%
  mutate(
    Entrez_ID = as.character(Entrez_ID),
    Concentration = case_when(
      str_detect(Condition, "0.1") ~ "0.1 µM",
      str_detect(Condition, "0.5") ~ "0.5 µM",
      TRUE ~ NA_character_
    )
  )

deg_mean <- deg_all %>%
  group_by(Entrez_ID, Concentration) %>%
  summarise(logFC = mean(logFC, na.rm = TRUE), .groups = "drop")

# ----------------- Merge with CorMotif TOP2B Genes -----------------
top2b_logfc <- inner_join(df_combined, deg_mean, by = c("Entrez_ID", "Concentration")) %>%
  mutate(
    Gene = mapIds(org.Hs.eg.db, keys = Entrez_ID,
                  column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"),
    Group = paste(Response_Group, Concentration)
  )

# ----------------- Define order -----------------
group_levels <- c(
  "Non response",
  "CX-DOX mid-late response",
  "DOX-specific response",
  "DOX only mid-late response",
  "CX total + DOX early response",
  "DOX early + CX-DOX mid-late response"
)

top2b_logfc$Response_Group <- factor(top2b_logfc$Response_Group, levels = group_levels)

# ----------------- Wilcoxon Test -----------------
comparison_map <- list(
  "CX-DOX mid-late response 0.1 µM" = "Non response 0.1 µM",
  "DOX only mid-late response 0.1 µM" = "Non response 0.1 µM",
  "DOX-specific response 0.5 µM" = "Non response 0.5 µM",
  "DOX only mid-late response 0.5 µM" = "Non response 0.5 µM",
  "CX total + DOX early response 0.5 µM" = "Non response 0.5 µM",
  "DOX early + CX-DOX mid-late response 0.5 µM" = "Non response 0.5 µM"
)

star_df <- lapply(names(comparison_map), function(resp_group) {
  control_group <- comparison_map[[resp_group]]
  resp_vals <- top2b_logfc$logFC[top2b_logfc$Group == resp_group]
  control_vals <- top2b_logfc$logFC[top2b_logfc$Group == control_group]

  if (length(resp_vals) >= 2 && length(control_vals) >= 2) {
    test_result <- wilcox.test(resp_vals, control_vals)
    pval <- test_result$p.value
    if (pval < 0.05) {
      label <- case_when(
        pval < 0.001 ~ "***",
        pval < 0.01 ~ "**",
        TRUE ~ "*"
      )
      y_pos <- max(resp_vals, na.rm = TRUE) + 0.5
      return(data.frame(Group = resp_group, y_pos = y_pos, label = label, P_Value = signif(pval, 4)))
    }
  }
  return(NULL)
}) %>% bind_rows()

# ----------------- Prepare Annotation -----------------
label_data <- star_df %>%
  separate(Group, into = c("Response_Group", "Concentration"), sep = " (?=0\\.)") %>%
  mutate(Response_Group = factor(Response_Group, levels = group_levels))

# ----------------- Plot -----------------
ggplot(top2b_logfc, aes(x = Response_Group, y = logFC, fill = Response_Group)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7) +
  geom_text(data = label_data, aes(x = Response_Group, y = y_pos, label = label),
            inherit.aes = FALSE, size = 5, color = "black") +
  facet_wrap(~Concentration, scales = "free_x") +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Mean logFC of TOP2B Target Genes\nAcross CorMotif Response Groups",
    x = "Response Group",
    y = "Mean Log Fold Change",
    fill = "Response Group"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 13),
    strip.text = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)
  )






############ TOP2A TARGET GENE ANALYSIS (ChIP Atlas) #######################

# ----------------- Load Libraries -----------------
library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)

# ----------------- Load TOP2A Entrez IDs -----------------
top2a_entrez_ids <- read_csv("data/TOP2A_target_mapped.csv", show_col_types = FALSE) %>%
  pull(Entrez_ID) %>%
  unique() %>%
  as.character()

# ----------------- Load CorMotif Groups -----------------
# 0.1 µM
prob_1_0.1 <- as.character(read.csv("data/prob_1_0.1.csv")$Entrez_ID)
prob_2_0.1 <- as.character(read.csv("data/prob_2_0.1.csv")$Entrez_ID)
prob_3_0.1 <- as.character(read.csv("data/prob_3_0.1.csv")$Entrez_ID)

# 0.5 µM
prob_1_0.5 <- as.character(read.csv("data/prob_1_0.5.csv")$Entrez_ID)
prob_2_0.5 <- as.character(read.csv("data/prob_2_0.5.csv")$Entrez_ID)
prob_3_0.5 <- as.character(read.csv("data/prob_3_0.5.csv")$Entrez_ID)
prob_4_0.5 <- as.character(read.csv("data/prob_4_0.5.csv")$Entrez_ID)
prob_5_0.5 <- as.character(read.csv("data/prob_5_0.5.csv")$Entrez_ID)

# ----------------- Annotate Response Groups -----------------
df_0.1 <- data.frame(Entrez_ID = unique(c(prob_1_0.1, prob_2_0.1, prob_3_0.1))) %>%
  filter(Entrez_ID %in% top2a_entrez_ids) %>%
  mutate(
    Response_Group = case_when(
      Entrez_ID %in% prob_1_0.1 ~ "Non response",
      Entrez_ID %in% prob_2_0.1 ~ "CX-DOX mid-late response",
      Entrez_ID %in% prob_3_0.1 ~ "DOX only mid-late response"
    ),
    Concentration = "0.1 µM"
  )

df_0.5 <- data.frame(Entrez_ID = unique(c(prob_1_0.5, prob_2_0.5, prob_3_0.5, prob_4_0.5, prob_5_0.5))) %>%
  filter(Entrez_ID %in% top2a_entrez_ids) %>%
  mutate(
    Response_Group = case_when(
      Entrez_ID %in% prob_1_0.5 ~ "Non response",
      Entrez_ID %in% prob_2_0.5 ~ "DOX-specific response",
      Entrez_ID %in% prob_3_0.5 ~ "DOX only mid-late response",
      Entrez_ID %in% prob_4_0.5 ~ "CX total + DOX early response",
      Entrez_ID %in% prob_5_0.5 ~ "DOX early + CX-DOX mid-late response"
    ),
    Concentration = "0.5 µM"
  )

df_combined <- bind_rows(df_0.1, df_0.5)

# ----------------- Load and Aggregate DEG Data -----------------
deg_files <- list(
  "CX_0.1_3" = "data/DEGs/Toptable_CX_0.1_3.csv",
  "CX_0.1_24" = "data/DEGs/Toptable_CX_0.1_24.csv",
  "CX_0.1_48" = "data/DEGs/Toptable_CX_0.1_48.csv",
  "CX_0.5_3" = "data/DEGs/Toptable_CX_0.5_3.csv",
  "CX_0.5_24" = "data/DEGs/Toptable_CX_0.5_24.csv",
  "CX_0.5_48" = "data/DEGs/Toptable_CX_0.5_48.csv",
  "DOX_0.1_3" = "data/DEGs/Toptable_DOX_0.1_3.csv",
  "DOX_0.1_24" = "data/DEGs/Toptable_DOX_0.1_24.csv",
  "DOX_0.1_48" = "data/DEGs/Toptable_DOX_0.1_48.csv",
  "DOX_0.5_3" = "data/DEGs/Toptable_DOX_0.5_3.csv",
  "DOX_0.5_24" = "data/DEGs/Toptable_DOX_0.5_24.csv",
  "DOX_0.5_48" = "data/DEGs/Toptable_DOX_0.5_48.csv"
)

deg_all <- map2_dfr(deg_files, names(deg_files), ~{
  read_csv(.x, show_col_types = FALSE) %>%
    mutate(Condition = .y)
}) %>%
  mutate(
    Entrez_ID = as.character(Entrez_ID),
    Concentration = case_when(
      str_detect(Condition, "0.1") ~ "0.1 µM",
      str_detect(Condition, "0.5") ~ "0.5 µM",
      TRUE ~ NA_character_
    )
  )

deg_mean <- deg_all %>%
  group_by(Entrez_ID, Concentration) %>%
  summarise(logFC = mean(logFC, na.rm = TRUE), .groups = "drop")

# ----------------- Merge with CorMotif TOP2A Genes -----------------
top2a_logfc <- inner_join(df_combined, deg_mean, by = c("Entrez_ID", "Concentration")) %>%
  mutate(
    Gene = mapIds(org.Hs.eg.db, keys = Entrez_ID,
                  column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"),
    Group = paste(Response_Group, Concentration)
  )

# ----------------- Define Response Group Order -----------------
group_levels <- c(
  "Non response",
  "CX-DOX mid-late response",
  "DOX-specific response",
  "DOX only mid-late response",
  "CX total + DOX early response",
  "DOX early + CX-DOX mid-late response"
)

top2a_logfc$Response_Group <- factor(top2a_logfc$Response_Group, levels = group_levels)

# ----------------- Wilcoxon Test with Output -----------------
comparison_map <- list(
  "CX-DOX mid-late response 0.1 µM" = "Non response 0.1 µM",
  "DOX only mid-late response 0.1 µM" = "Non response 0.1 µM",
  "DOX-specific response 0.5 µM" = "Non response 0.5 µM",
  "DOX only mid-late response 0.5 µM" = "Non response 0.5 µM",
  "CX total + DOX early response 0.5 µM" = "Non response 0.5 µM",
  "DOX early + CX-DOX mid-late response 0.5 µM" = "Non response 0.5 µM"
)

wilcoxon_results <- lapply(names(comparison_map), function(resp_group) {
  control_group <- comparison_map[[resp_group]]
  resp_vals <- top2a_logfc$logFC[top2a_logfc$Group == resp_group]
  control_vals <- top2a_logfc$logFC[top2a_logfc$Group == control_group]

  if (length(resp_vals) >= 2 && length(control_vals) >= 2) {
    test_result <- wilcox.test(resp_vals, control_vals)
    pval <- test_result$p.value
    stat <- test_result$statistic
    label <- case_when(
      pval < 0.001 ~ "***",
      pval < 0.01  ~ "**",
      pval < 0.05  ~ "*",
      TRUE ~ ""
    )
    y_pos <- max(resp_vals, na.rm = TRUE) + 0.5
    return(data.frame(
      Response_Group = resp_group,
      Control_Group = control_group,
      W_statistic = as.numeric(stat),
      P_value = signif(pval, 4),
      Significance = label,
      y_pos = y_pos
    ))
  } else {
    return(data.frame(
      Response_Group = resp_group,
      Control_Group = control_group,
      W_statistic = NA,
      P_value = NA,
      Significance = "",
      y_pos = NA
    ))
  }
}) %>% bind_rows()

# ----------------- Prepare Annotation -----------------
label_data <- wilcoxon_results %>%
  separate(Response_Group, into = c("Response_Group", "Concentration"), sep = " (?=0\\.)") %>%
  mutate(Response_Group = factor(Response_Group, levels = group_levels))

# ----------------- Plot -----------------
ggplot(top2a_logfc, aes(x = Response_Group, y = logFC, fill = Response_Group)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7) +
  geom_text(data = label_data, aes(x = Response_Group, y = y_pos, label = Significance),
            inherit.aes = FALSE, size = 5, color = "black") +
  facet_wrap(~Concentration, scales = "free_x") +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Mean logFC of TOP2A Target Genes\nAcross CorMotif Response Groups (ChIP atlas)",
    x = "Response Group",
    y = "Mean Log Fold Change",
    fill = "Response Group"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 13),
    strip.text = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)
  )

# Optional: View or save results
print(wilcoxon_results)
# write_csv(wilcoxon_results, "TOP2A_Wilcoxon_Results.csv")



############ TOP2A LITERATURE TARGET GENE ANALYSIS (RPE-1 Cell Line) #######################

# ----------------- Load Libraries -----------------
library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)

# ----------------- Load TOP2A (RPE-1) Entrez IDs -----------------
top2a_lit_ids <- read_csv("data/TOP2A_target_lit_mapped.csv", show_col_types = FALSE) %>%
  pull(Entrez_ID) %>%
  unique() %>%
  as.character()

# ----------------- Load CorMotif Groups -----------------
# 0.1 µM
prob_1_0.1 <- as.character(read.csv("data/prob_1_0.1.csv")$Entrez_ID)
prob_2_0.1 <- as.character(read.csv("data/prob_2_0.1.csv")$Entrez_ID)
prob_3_0.1 <- as.character(read.csv("data/prob_3_0.1.csv")$Entrez_ID)

# 0.5 µM
prob_1_0.5 <- as.character(read.csv("data/prob_1_0.5.csv")$Entrez_ID)
prob_2_0.5 <- as.character(read.csv("data/prob_2_0.5.csv")$Entrez_ID)
prob_3_0.5 <- as.character(read.csv("data/prob_3_0.5.csv")$Entrez_ID)
prob_4_0.5 <- as.character(read.csv("data/prob_4_0.5.csv")$Entrez_ID)
prob_5_0.5 <- as.character(read.csv("data/prob_5_0.5.csv")$Entrez_ID)

# ----------------- Annotate Response Groups -----------------
df_0.1 <- data.frame(Entrez_ID = unique(c(prob_1_0.1, prob_2_0.1, prob_3_0.1))) %>%
  filter(Entrez_ID %in% top2a_lit_ids) %>%
  mutate(
    Response_Group = case_when(
      Entrez_ID %in% prob_1_0.1 ~ "Non response",
      Entrez_ID %in% prob_2_0.1 ~ "CX-DOX mid-late response",
      Entrez_ID %in% prob_3_0.1 ~ "DOX only mid-late response"
    ),
    Concentration = "0.1 µM"
  )

df_0.5 <- data.frame(Entrez_ID = unique(c(prob_1_0.5, prob_2_0.5, prob_3_0.5, prob_4_0.5, prob_5_0.5))) %>%
  filter(Entrez_ID %in% top2a_lit_ids) %>%
  mutate(
    Response_Group = case_when(
      Entrez_ID %in% prob_1_0.5 ~ "Non response",
      Entrez_ID %in% prob_2_0.5 ~ "DOX-specific response",
      Entrez_ID %in% prob_3_0.5 ~ "DOX only mid-late response",
      Entrez_ID %in% prob_4_0.5 ~ "CX total + DOX early response",
      Entrez_ID %in% prob_5_0.5 ~ "DOX early + CX-DOX mid-late response"
    ),
    Concentration = "0.5 µM"
  )

df_combined <- bind_rows(df_0.1, df_0.5)

# ----------------- Load and Aggregate DEG Data -----------------
deg_files <- list(
  "CX_0.1_3" = "data/DEGs/Toptable_CX_0.1_3.csv",
  "CX_0.1_24" = "data/DEGs/Toptable_CX_0.1_24.csv",
  "CX_0.1_48" = "data/DEGs/Toptable_CX_0.1_48.csv",
  "CX_0.5_3" = "data/DEGs/Toptable_CX_0.5_3.csv",
  "CX_0.5_24" = "data/DEGs/Toptable_CX_0.5_24.csv",
  "CX_0.5_48" = "data/DEGs/Toptable_CX_0.5_48.csv",
  "DOX_0.1_3" = "data/DEGs/Toptable_DOX_0.1_3.csv",
  "DOX_0.1_24" = "data/DEGs/Toptable_DOX_0.1_24.csv",
  "DOX_0.1_48" = "data/DEGs/Toptable_DOX_0.1_48.csv",
  "DOX_0.5_3" = "data/DEGs/Toptable_DOX_0.5_3.csv",
  "DOX_0.5_24" = "data/DEGs/Toptable_DOX_0.5_24.csv",
  "DOX_0.5_48" = "data/DEGs/Toptable_DOX_0.5_48.csv"
)

deg_all <- map2_dfr(deg_files, names(deg_files), ~{
  read_csv(.x, show_col_types = FALSE) %>%
    mutate(Condition = .y)
}) %>%
  mutate(
    Entrez_ID = as.character(Entrez_ID),
    Concentration = case_when(
      str_detect(Condition, "0.1") ~ "0.1 µM",
      str_detect(Condition, "0.5") ~ "0.5 µM",
      TRUE ~ NA_character_
    )
  )

deg_mean <- deg_all %>%
  group_by(Entrez_ID, Concentration) %>%
  summarise(logFC = mean(logFC, na.rm = TRUE), .groups = "drop")

# ----------------- Merge with CorMotif Genes -----------------
top2a_lit_logfc <- inner_join(df_combined, deg_mean, by = c("Entrez_ID", "Concentration")) %>%
  mutate(
    Gene = mapIds(org.Hs.eg.db, keys = Entrez_ID,
                  column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"),
    Group = paste(Response_Group, Concentration)
  )

# ----------------- Define Response Group Order -----------------
group_levels <- c(
  "Non response",
  "CX-DOX mid-late response",
  "DOX-specific response",
  "DOX only mid-late response",
  "CX total + DOX early response",
  "DOX early + CX-DOX mid-late response"
)

top2a_lit_logfc$Response_Group <- factor(top2a_lit_logfc$Response_Group, levels = group_levels)

# ----------------- Wilcoxon Test with Output -----------------
comparison_map <- list(
  "CX-DOX mid-late response 0.1 µM" = "Non response 0.1 µM",
  "DOX only mid-late response 0.1 µM" = "Non response 0.1 µM",
  "DOX-specific response 0.5 µM" = "Non response 0.5 µM",
  "DOX only mid-late response 0.5 µM" = "Non response 0.5 µM",
  "CX total + DOX early response 0.5 µM" = "Non response 0.5 µM",
  "DOX early + CX-DOX mid-late response 0.5 µM" = "Non response 0.5 µM"
)

wilcoxon_results <- lapply(names(comparison_map), function(resp_group) {
  control_group <- comparison_map[[resp_group]]
  resp_vals <- top2a_lit_logfc$logFC[top2a_lit_logfc$Group == resp_group]
  control_vals <- top2a_lit_logfc$logFC[top2a_lit_logfc$Group == control_group]

  if (length(resp_vals) >= 2 && length(control_vals) >= 2) {
    test_result <- wilcox.test(resp_vals, control_vals)
    pval <- test_result$p.value
    stat <- test_result$statistic
    label <- case_when(
      pval < 0.001 ~ "***",
      pval < 0.01  ~ "**",
      pval < 0.05  ~ "*",
      TRUE ~ ""
    )
    y_pos <- max(resp_vals, na.rm = TRUE) + 0.5
    return(data.frame(
      Response_Group = resp_group,
      Control_Group = control_group,
      W_statistic = as.numeric(stat),
      P_value = signif(pval, 4),
      Significance = label,
      y_pos = y_pos
    ))
  } else {
    return(data.frame(
      Response_Group = resp_group,
      Control_Group = control_group,
      W_statistic = NA,
      P_value = NA,
      Significance = "",
      y_pos = NA
    ))
  }
}) %>% bind_rows()

# ----------------- Prepare Annotation -----------------
label_data <- wilcoxon_results %>%
  separate(Response_Group, into = c("Response_Group", "Concentration"), sep = " (?=0\\.)") %>%
  mutate(Response_Group = factor(Response_Group, levels = group_levels))

# ----------------- Plot -----------------
ggplot(top2a_lit_logfc, aes(x = Response_Group, y = logFC, fill = Response_Group)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7) +
  geom_text(data = label_data, aes(x = Response_Group, y = y_pos, label = Significance),
            inherit.aes = FALSE, size = 5, color = "black") +
  facet_wrap(~Concentration, scales = "free_x") +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Mean logFC of Literature-Based TOP2A Targets (RPE-1)\nAcross CorMotif Response Groups",
    x = "Response Group",
    y = "Mean Log Fold Change",
    fill = "Response Group"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 13),
    strip.text = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)
  )

# Optional: View or save results
print(wilcoxon_results)
# write_csv(wilcoxon_results, "TOP2A_Lit_Wilcoxon_Results.csv")




############ AC CARDIOTOXICITY GENE ANALYSIS (CorMotif-based) #######################

# ----------------- Load Libraries -----------------
library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)

# ----------------- AC Cardiotoxicity Entrez IDs -----------------
entrez_ids <- c(
  6272, 8029, 11128, 79899, 54477, 121665, 5095, 22863, 57161, 4692,
  8214, 23151, 56606, 108, 22999, 56895, 9603, 3181, 4023, 10499,
  92949, 4363, 10057, 5243, 5244, 5880, 1535, 2950, 847, 5447,
  3038, 3077, 4846, 3958, 23327, 29899, 23155, 80856, 55020, 78996,
  23262, 150383, 9620, 79730, 344595, 5066, 6251, 3482, 9588, 339416,
  7292, 55157, 87769, 23409, 720, 3107, 54535, 1590, 80059, 7991,
  57110, 8803, 323, 54826, 5916, 23371, 283337, 64078, 80010, 1933,
  10818, 51020
) %>% as.character()

# ----------------- Load CorMotif Groups -----------------
# 0.1 µM
prob_1_0.1 <- as.character(read.csv("data/prob_1_0.1.csv")$Entrez_ID)
prob_2_0.1 <- as.character(read.csv("data/prob_2_0.1.csv")$Entrez_ID)
prob_3_0.1 <- as.character(read.csv("data/prob_3_0.1.csv")$Entrez_ID)

# 0.5 µM
prob_1_0.5 <- as.character(read.csv("data/prob_1_0.5.csv")$Entrez_ID)
prob_2_0.5 <- as.character(read.csv("data/prob_2_0.5.csv")$Entrez_ID)
prob_3_0.5 <- as.character(read.csv("data/prob_3_0.5.csv")$Entrez_ID)
prob_4_0.5 <- as.character(read.csv("data/prob_4_0.5.csv")$Entrez_ID)
prob_5_0.5 <- as.character(read.csv("data/prob_5_0.5.csv")$Entrez_ID)

# ----------------- Annotate Response Groups -----------------
df_0.1 <- data.frame(Entrez_ID = unique(c(prob_1_0.1, prob_2_0.1, prob_3_0.1))) %>%
  filter(Entrez_ID %in% entrez_ids) %>%
  mutate(
    Response_Group = case_when(
      Entrez_ID %in% prob_1_0.1 ~ "Non response",
      Entrez_ID %in% prob_2_0.1 ~ "CX-DOX mid-late response",
      Entrez_ID %in% prob_3_0.1 ~ "DOX only mid-late response"
    ),
    Concentration = "0.1 µM"
  )

df_0.5 <- data.frame(Entrez_ID = unique(c(prob_1_0.5, prob_2_0.5, prob_3_0.5, prob_4_0.5, prob_5_0.5))) %>%
  filter(Entrez_ID %in% entrez_ids) %>%
  mutate(
    Response_Group = case_when(
      Entrez_ID %in% prob_1_0.5 ~ "Non response",
      Entrez_ID %in% prob_2_0.5 ~ "DOX-specific response",
      Entrez_ID %in% prob_3_0.5 ~ "DOX only mid-late response",
      Entrez_ID %in% prob_4_0.5 ~ "CX total + DOX early response",
      Entrez_ID %in% prob_5_0.5 ~ "DOX early + CX-DOX mid-late response"
    ),
    Concentration = "0.5 µM"
  )

df_combined <- bind_rows(df_0.1, df_0.5)

# ----------------- Load and Aggregate DEG Data -----------------
deg_files <- list(
  "CX_0.1_3" = "data/DEGs/Toptable_CX_0.1_3.csv",
  "CX_0.1_24" = "data/DEGs/Toptable_CX_0.1_24.csv",
  "CX_0.1_48" = "data/DEGs/Toptable_CX_0.1_48.csv",
  "CX_0.5_3" = "data/DEGs/Toptable_CX_0.5_3.csv",
  "CX_0.5_24" = "data/DEGs/Toptable_CX_0.5_24.csv",
  "CX_0.5_48" = "data/DEGs/Toptable_CX_0.5_48.csv",
  "DOX_0.1_3" = "data/DEGs/Toptable_DOX_0.1_3.csv",
  "DOX_0.1_24" = "data/DEGs/Toptable_DOX_0.1_24.csv",
  "DOX_0.1_48" = "data/DEGs/Toptable_DOX_0.1_48.csv",
  "DOX_0.5_3" = "data/DEGs/Toptable_DOX_0.5_3.csv",
  "DOX_0.5_24" = "data/DEGs/Toptable_DOX_0.5_24.csv",
  "DOX_0.5_48" = "data/DEGs/Toptable_DOX_0.5_48.csv"
)

deg_all <- map2_dfr(deg_files, names(deg_files), ~{
  read_csv(.x, show_col_types = FALSE) %>%
    mutate(Condition = .y)
}) %>%
  mutate(
    Entrez_ID = as.character(Entrez_ID),
    Concentration = case_when(
      str_detect(Condition, "0.1") ~ "0.1 µM",
      str_detect(Condition, "0.5") ~ "0.5 µM",
      TRUE ~ NA_character_
    )
  )

deg_mean <- deg_all %>%
  group_by(Entrez_ID, Concentration) %>%
  summarise(logFC = mean(logFC, na.rm = TRUE), .groups = "drop")

# ----------------- Merge with CorMotif Genes -----------------
ac_logfc <- inner_join(df_combined, deg_mean, by = c("Entrez_ID", "Concentration")) %>%
  mutate(
    Gene = mapIds(org.Hs.eg.db, keys = Entrez_ID,
                  column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"),
    Group = paste(Response_Group, Concentration)
  )

# ----------------- Define Group Order -----------------
group_levels <- c(
  "Non response",
  "CX-DOX mid-late response",
  "DOX-specific response",
  "DOX only mid-late response",
  "CX total + DOX early response",
  "DOX early + CX-DOX mid-late response"
)

ac_logfc$Response_Group <- factor(ac_logfc$Response_Group, levels = group_levels)

# ----------------- Wilcoxon Test -----------------
comparison_map <- list(
  "CX-DOX mid-late response 0.1 µM" = "Non response 0.1 µM",
  "DOX only mid-late response 0.1 µM" = "Non response 0.1 µM",
  "DOX-specific response 0.5 µM" = "Non response 0.5 µM",
  "DOX only mid-late response 0.5 µM" = "Non response 0.5 µM",
  "CX total + DOX early response 0.5 µM" = "Non response 0.5 µM",
  "DOX early + CX-DOX mid-late response 0.5 µM" = "Non response 0.5 µM"
)

wilcoxon_results <- lapply(names(comparison_map), function(resp_group) {
  control_group <- comparison_map[[resp_group]]
  resp_vals <- ac_logfc$logFC[ac_logfc$Group == resp_group]
  control_vals <- ac_logfc$logFC[ac_logfc$Group == control_group]

  if (length(resp_vals) >= 2 && length(control_vals) >= 2) {
    test_result <- wilcox.test(resp_vals, control_vals)
    pval <- test_result$p.value
    stat <- test_result$statistic
    label <- case_when(
      pval < 0.001 ~ "***",
      pval < 0.01  ~ "**",
      pval < 0.05  ~ "*",
      TRUE ~ ""
    )
    y_pos <- max(resp_vals, na.rm = TRUE) + 0.5
    return(data.frame(
      Group = resp_group,
      y_pos = y_pos,
      label = label,
      P_Value = signif(pval, 4)
    ))
  } else {
    return(NULL)
  }
}) %>% bind_rows()

# ----------------- Prepare Annotation -----------------
label_data <- wilcoxon_results %>%
  separate(Group, into = c("Response_Group", "Concentration"), sep = " (?=0\\.)") %>%
  mutate(Response_Group = factor(Response_Group, levels = group_levels))

# ----------------- Plot -----------------
ggplot(ac_logfc, aes(x = Response_Group, y = logFC, fill = Response_Group)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7) +
  geom_text(data = label_data, aes(x = Response_Group, y = y_pos, label = label),
            inherit.aes = FALSE, size = 5, color = "black") +
  facet_wrap(~Concentration, scales = "free_x") +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Mean logFC of AC Cardiotoxicity Genes\nAcross CorMotif Response Groups",
    x = "Response Group",
    y = "Mean Log Fold Change",
    fill = "Response Group"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 13),
    strip.text = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)
  )



############ DOX CARDIOTOXICITY GENE ANALYSIS #######################

# ----------------- Load Libraries -----------------
library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)

# ----------------- DOX Cardiotoxicity Entrez IDs -----------------
entrez_ids <- c(
  847, 873, 2064, 2878, 2944, 3038, 4846, 51196, 5880, 6687,
  7799, 4292, 5916, 3077, 51310, 9154, 64078, 5244, 10057, 10060,
  89845, 56853, 4625, 1573, 79890
) %>% as.character()

# ----------------- Load CorMotif Groups -----------------
# 0.1 µM
prob_1_0.1 <- as.character(read.csv("data/prob_1_0.1.csv")$Entrez_ID)
prob_2_0.1 <- as.character(read.csv("data/prob_2_0.1.csv")$Entrez_ID)
prob_3_0.1 <- as.character(read.csv("data/prob_3_0.1.csv")$Entrez_ID)

# 0.5 µM
prob_1_0.5 <- as.character(read.csv("data/prob_1_0.5.csv")$Entrez_ID)
prob_2_0.5 <- as.character(read.csv("data/prob_2_0.5.csv")$Entrez_ID)
prob_3_0.5 <- as.character(read.csv("data/prob_3_0.5.csv")$Entrez_ID)
prob_4_0.5 <- as.character(read.csv("data/prob_4_0.5.csv")$Entrez_ID)
prob_5_0.5 <- as.character(read.csv("data/prob_5_0.5.csv")$Entrez_ID)

# ----------------- Annotate Response Groups -----------------
df_0.1 <- data.frame(Entrez_ID = unique(c(prob_1_0.1, prob_2_0.1, prob_3_0.1))) %>%
  filter(Entrez_ID %in% entrez_ids) %>%
  mutate(
    Response_Group = case_when(
      Entrez_ID %in% prob_1_0.1 ~ "Non response",
      Entrez_ID %in% prob_2_0.1 ~ "CX-DOX mid-late response",
      Entrez_ID %in% prob_3_0.1 ~ "DOX only mid-late response"
    ),
    Concentration = "0.1 µM"
  )

df_0.5 <- data.frame(Entrez_ID = unique(c(prob_1_0.5, prob_2_0.5, prob_3_0.5, prob_4_0.5, prob_5_0.5))) %>%
  filter(Entrez_ID %in% entrez_ids) %>%
  mutate(
    Response_Group = case_when(
      Entrez_ID %in% prob_1_0.5 ~ "Non response",
      Entrez_ID %in% prob_2_0.5 ~ "DOX-specific response",
      Entrez_ID %in% prob_3_0.5 ~ "DOX only mid-late response",
      Entrez_ID %in% prob_4_0.5 ~ "CX total + DOX early response",
      Entrez_ID %in% prob_5_0.5 ~ "DOX early + CX-DOX mid-late response"
    ),
    Concentration = "0.5 µM"
  )

df_combined <- bind_rows(df_0.1, df_0.5)

# ----------------- Load and Aggregate DEG Data -----------------
deg_files <- list(
  "CX_0.1_3" = "data/DEGs/Toptable_CX_0.1_3.csv",
  "CX_0.1_24" = "data/DEGs/Toptable_CX_0.1_24.csv",
  "CX_0.1_48" = "data/DEGs/Toptable_CX_0.1_48.csv",
  "CX_0.5_3" = "data/DEGs/Toptable_CX_0.5_3.csv",
  "CX_0.5_24" = "data/DEGs/Toptable_CX_0.5_24.csv",
  "CX_0.5_48" = "data/DEGs/Toptable_CX_0.5_48.csv",
  "DOX_0.1_3" = "data/DEGs/Toptable_DOX_0.1_3.csv",
  "DOX_0.1_24" = "data/DEGs/Toptable_DOX_0.1_24.csv",
  "DOX_0.1_48" = "data/DEGs/Toptable_DOX_0.1_48.csv",
  "DOX_0.5_3" = "data/DEGs/Toptable_DOX_0.5_3.csv",
  "DOX_0.5_24" = "data/DEGs/Toptable_DOX_0.5_24.csv",
  "DOX_0.5_48" = "data/DEGs/Toptable_DOX_0.5_48.csv"
)

deg_all <- map2_dfr(deg_files, names(deg_files), ~{
  read_csv(.x, show_col_types = FALSE) %>%
    mutate(Condition = .y)
}) %>%
  mutate(
    Entrez_ID = as.character(Entrez_ID),
    Concentration = case_when(
      str_detect(Condition, "0.1") ~ "0.1 µM",
      str_detect(Condition, "0.5") ~ "0.5 µM",
      TRUE ~ NA_character_
    )
  )

deg_mean <- deg_all %>%
  group_by(Entrez_ID, Concentration) %>%
  summarise(logFC = mean(logFC, na.rm = TRUE), .groups = "drop")

# ----------------- Merge with CorMotif Genes -----------------
dox_cardio_logfc <- inner_join(df_combined, deg_mean, by = c("Entrez_ID", "Concentration")) %>%
  mutate(
    Gene = mapIds(org.Hs.eg.db, keys = Entrez_ID,
                  column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"),
    Group = paste(Response_Group, Concentration)
  )

# ----------------- Define Response Group Order -----------------
group_levels <- c(
  "Non response",
  "CX-DOX mid-late response",
  "DOX-specific response",
  "DOX only mid-late response",
  "CX total + DOX early response",
  "DOX early + CX-DOX mid-late response"
)

dox_cardio_logfc$Response_Group <- factor(dox_cardio_logfc$Response_Group, levels = group_levels)

# ----------------- Wilcoxon Test -----------------
comparison_map <- list(
  "CX-DOX mid-late response 0.1 µM" = "Non response 0.1 µM",
  "DOX only mid-late response 0.1 µM" = "Non response 0.1 µM",
  "DOX-specific response 0.5 µM" = "Non response 0.5 µM",
  "DOX only mid-late response 0.5 µM" = "Non response 0.5 µM",
  "CX total + DOX early response 0.5 µM" = "Non response 0.5 µM",
  "DOX early + CX-DOX mid-late response 0.5 µM" = "Non response 0.5 µM"
)

wilcoxon_results <- lapply(names(comparison_map), function(resp_group) {
  control_group <- comparison_map[[resp_group]]
  resp_vals <- dox_cardio_logfc$logFC[dox_cardio_logfc$Group == resp_group]
  control_vals <- dox_cardio_logfc$logFC[dox_cardio_logfc$Group == control_group]

  if (length(resp_vals) >= 2 && length(control_vals) >= 2) {
    test_result <- wilcox.test(resp_vals, control_vals)
    pval <- test_result$p.value
    stat <- test_result$statistic
    label <- case_when(
      pval < 0.001 ~ "***",
      pval < 0.01  ~ "**",
      pval < 0.05  ~ "*",
      TRUE ~ ""
    )
    y_pos <- max(resp_vals, na.rm = TRUE) + 0.5
    return(data.frame(
      Group = resp_group,
      y_pos = y_pos,
      label = label,
      P_Value = signif(pval, 4)
    ))
  } else {
    return(NULL)
  }
}) %>% bind_rows()

# ----------------- Prepare Annotation -----------------
label_data <- wilcoxon_results %>%
  separate(Group, into = c("Response_Group", "Concentration"), sep = " (?=0\\.)") %>%
  mutate(Response_Group = factor(Response_Group, levels = group_levels))

# ----------------- Plot -----------------
ggplot(dox_cardio_logfc, aes(x = Response_Group, y = logFC, fill = Response_Group)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7) +
  geom_text(data = label_data, aes(x = Response_Group, y = y_pos, label = label),
            inherit.aes = FALSE, size = 5, color = "black") +
  facet_wrap(~Concentration, scales = "free_x") +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Mean logFC of DOX Cardiotoxicity Genes\nAcross CorMotif Response Groups",
    x = "Response Group",
    y = "Mean Log Fold Change",
    fill = "Response Group"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 13),
    strip.text = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)
  )







############ HEART-SPECIFIC GENE ANALYSIS #######################

# ----------------- Load Libraries -----------------
library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)

# ----------------- Load Heart-Specific Genes -----------------
heart_genes <- read.csv("data/Human_Heart_Genes.csv", stringsAsFactors = FALSE)
heart_genes$Entrez_ID <- mapIds(
  org.Hs.eg.db,
  keys = heart_genes$Gene,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)
entrez_ids <- na.omit(heart_genes$Entrez_ID) %>% as.character()

# ----------------- Load CorMotif Groups -----------------
# 0.1 µM
prob_1_0.1 <- as.character(read.csv("data/prob_1_0.1.csv")$Entrez_ID)
prob_2_0.1 <- as.character(read.csv("data/prob_2_0.1.csv")$Entrez_ID)
prob_3_0.1 <- as.character(read.csv("data/prob_3_0.1.csv")$Entrez_ID)

# 0.5 µM
prob_1_0.5 <- as.character(read.csv("data/prob_1_0.5.csv")$Entrez_ID)
prob_2_0.5 <- as.character(read.csv("data/prob_2_0.5.csv")$Entrez_ID)
prob_3_0.5 <- as.character(read.csv("data/prob_3_0.5.csv")$Entrez_ID)
prob_4_0.5 <- as.character(read.csv("data/prob_4_0.5.csv")$Entrez_ID)
prob_5_0.5 <- as.character(read.csv("data/prob_5_0.5.csv")$Entrez_ID)

# ----------------- Annotate Response Groups -----------------
df_0.1 <- data.frame(Entrez_ID = unique(c(prob_1_0.1, prob_2_0.1, prob_3_0.1))) %>%
  filter(Entrez_ID %in% entrez_ids) %>%
  mutate(
    Response_Group = case_when(
      Entrez_ID %in% prob_1_0.1 ~ "Non response",
      Entrez_ID %in% prob_2_0.1 ~ "CX-DOX mid-late response",
      Entrez_ID %in% prob_3_0.1 ~ "DOX only mid-late response"
    ),
    Concentration = "0.1 µM"
  )

df_0.5 <- data.frame(Entrez_ID = unique(c(prob_1_0.5, prob_2_0.5, prob_3_0.5, prob_4_0.5, prob_5_0.5))) %>%
  filter(Entrez_ID %in% entrez_ids) %>%
  mutate(
    Response_Group = case_when(
      Entrez_ID %in% prob_1_0.5 ~ "Non response",
      Entrez_ID %in% prob_2_0.5 ~ "DOX-specific response",
      Entrez_ID %in% prob_3_0.5 ~ "DOX only mid-late response",
      Entrez_ID %in% prob_4_0.5 ~ "CX total + DOX early response",
      Entrez_ID %in% prob_5_0.5 ~ "DOX early + CX-DOX mid-late response"
    ),
    Concentration = "0.5 µM"
  )

df_combined <- bind_rows(df_0.1, df_0.5)

# ----------------- Load and Aggregate DEG Data -----------------
deg_files <- list(
  "CX_0.1_3" = "data/DEGs/Toptable_CX_0.1_3.csv",
  "CX_0.1_24" = "data/DEGs/Toptable_CX_0.1_24.csv",
  "CX_0.1_48" = "data/DEGs/Toptable_CX_0.1_48.csv",
  "CX_0.5_3" = "data/DEGs/Toptable_CX_0.5_3.csv",
  "CX_0.5_24" = "data/DEGs/Toptable_CX_0.5_24.csv",
  "CX_0.5_48" = "data/DEGs/Toptable_CX_0.5_48.csv",
  "DOX_0.1_3" = "data/DEGs/Toptable_DOX_0.1_3.csv",
  "DOX_0.1_24" = "data/DEGs/Toptable_DOX_0.1_24.csv",
  "DOX_0.1_48" = "data/DEGs/Toptable_DOX_0.1_48.csv",
  "DOX_0.5_3" = "data/DEGs/Toptable_DOX_0.5_3.csv",
  "DOX_0.5_24" = "data/DEGs/Toptable_DOX_0.5_24.csv",
  "DOX_0.5_48" = "data/DEGs/Toptable_DOX_0.5_48.csv"
)

deg_all <- map2_dfr(deg_files, names(deg_files), ~{
  read_csv(.x, show_col_types = FALSE) %>% mutate(Condition = .y)
}) %>%
  mutate(
    Entrez_ID = as.character(Entrez_ID),
    Concentration = case_when(
      str_detect(Condition, "0.1") ~ "0.1 µM",
      str_detect(Condition, "0.5") ~ "0.5 µM"
    )
  )

deg_mean <- deg_all %>%
  group_by(Entrez_ID, Concentration) %>%
  summarise(logFC = mean(logFC, na.rm = TRUE), .groups = "drop")

# ----------------- Merge with Group Annotations -----------------
heart_logfc <- inner_join(df_combined, deg_mean, by = c("Entrez_ID", "Concentration")) %>%
  mutate(
    Gene = mapIds(org.Hs.eg.db, keys = Entrez_ID,
                  column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"),
    Group = paste(Response_Group, Concentration)
  )

# ----------------- Define Response Group Order -----------------
group_levels <- c(
  "Non response",
  "CX-DOX mid-late response",
  "DOX-specific response",
  "DOX only mid-late response",
  "CX total + DOX early response",
  "DOX early + CX-DOX mid-late response"
)

heart_logfc$Response_Group <- factor(heart_logfc$Response_Group, levels = group_levels)

# ----------------- Wilcoxon Test -----------------
comparison_map <- list(
  "CX-DOX mid-late response 0.1 µM" = "Non response 0.1 µM",
  "DOX only mid-late response 0.1 µM" = "Non response 0.1 µM",
  "DOX-specific response 0.5 µM" = "Non response 0.5 µM",
  "DOX only mid-late response 0.5 µM" = "Non response 0.5 µM",
  "CX total + DOX early response 0.5 µM" = "Non response 0.5 µM",
  "DOX early + CX-DOX mid-late response 0.5 µM" = "Non response 0.5 µM"
)

wilcoxon_results <- lapply(names(comparison_map), function(resp_group) {
  control_group <- comparison_map[[resp_group]]
  resp_vals <- heart_logfc$logFC[heart_logfc$Group == resp_group]
  control_vals <- heart_logfc$logFC[heart_logfc$Group == control_group]

  if (length(resp_vals) >= 2 && length(control_vals) >= 2) {
    test_result <- wilcox.test(resp_vals, control_vals)
    pval <- test_result$p.value
    label <- case_when(
      pval < 0.001 ~ "***",
      pval < 0.01  ~ "**",
      pval < 0.05  ~ "*",
      TRUE ~ ""
    )
    y_pos <- max(resp_vals, na.rm = TRUE) + 0.5
    return(data.frame(
      Group = resp_group,
      y_pos = y_pos,
      label = label,
      P_Value = signif(pval, 4)
    ))
  }
  return(NULL)
}) %>% bind_rows()

# ----------------- Annotation Data -----------------
label_data <- wilcoxon_results %>%
  separate(Group, into = c("Response_Group", "Concentration"), sep = " (?=0\\.)") %>%
  mutate(Response_Group = factor(Response_Group, levels = group_levels))

# ----------------- Plot -----------------
ggplot(heart_logfc, aes(x = Response_Group, y = logFC, fill = Response_Group)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7) +
  geom_text(data = label_data, aes(x = Response_Group, y = y_pos, label = label),
            inherit.aes = FALSE, size = 5, color = "black") +
  facet_wrap(~Concentration, scales = "free_x") +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Mean logFC of Heart-Specific Genes\nAcross CorMotif Response Groups",
    x = "Response Group",
    y = "Mean Log Fold Change",
    fill = "Response Group"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 13),
    strip.text = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)
  )


############ ATRIAL FIBRILLATION GENE ANALYSIS #######################

# ----------------- Load Libraries -----------------
library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)

# ----------------- AF Entrez IDs -----------------
entrez_ids <- c(
  3425, 53834, 9464, 79658, 55876, 3781, 9472, 11017, 64395, 23353,
  57794, 79591, 2702, 51306, 857, 10418, 9644, 115509, 89797, 10728,
  282996, 8497, 5695, 9026, 57666, 143384, 5883, 90102, 84909, 221037,
  29119, 79933, 391, 81575, 8738, 5396, 26010, 744, 56981, 80724,
  5781, 6910, 4598, 79068, 2250, 56978, 3752, 5208, 9736, 26778,
  80315, 64374, 79159, 11078, 5308, 3782, 10654, 79776, 6310, 8082,
  2995, 23224, 79741, 5708, 26136, 9967, 84033, 9148, 29982, 81488,
  1788, 9208, 2908, 1021, 5339, 5727, 23118, 10090, 57585, 221154,
  1073, 29998, 817, 858, 3486, 5747, 2104, 55105, 4194, 3570,
  1607, 6586, 5570, 2070, 27145, 4646, 3159, 9748, 79991, 221035,
  5033, 80212, 196385, 79600, 11046, 54014, 151636, 64778, 92597, 3983,
  7090, 9267, 7273, 26018, 9948, 2131, 55870, 5523, 5318, 6239,
  3480, 55777, 2066, 11165, 2028, 57619, 2690, 23414, 55818, 84700,
  28965, 80204, 463, 5108, 222553, 387119, 28981, 11113, 10301, 995,
  23030, 2697, 996, 339500, 54805, 9960, 91404, 145781, 100820829, 84542,
  2176, 51684, 9513, 7473, 4666, 23150, 5915, 5062, 4016, 23039,
  159686, 1839, 5201, 93166, 64753, 29959, 5496, 23245, 5069, 56916,
  92344, 23092, 3992, 9415, 10554, 7456, 9570, 57178, 23143, 161176,
  5424, 2034, 10277, 11278, 79803, 6653, 4756, 132660, 5430, 9031,
  57158, 285761, 8110, 387700, 1829, 4126, 7323, 51308, 7332, 6598,
  3757, 6187, 6660, 10529, 6920, 115286, 8451, 8943, 4137, 7514,
  4801, 9709, 23177, 8671, 29915, 26207, 3680, 490, 493856, 58489,
  54897, 4625, 22955, 84952, 10221, 2263, 84641, 4892, 1026, 84650,
  8382, 221656, 2969, 144453, 117177, 2626, 8476, 161882, 51807, 4624,
  9612, 55795, 51043, 144348, 51232, 8729, 3899, 11155, 23316, 79006,
  146330, 6403, 6331, 4300, 427, 845, 3313, 113622, 5789, 376132,
  29841, 8462, 203859, 401397, 5506, 55521, 5819, 4059, 6934, 57727,
  23066, 79568, 83478, 10087, 9586, 222194, 10466, 10499, 58499, 79720,
  4772, 10794, 125919, 602, 27332, 59345, 340359, 3705, 64710, 57801,
  10818, 143684, 149281, 55013, 23095, 93649, 84034, 23347, 440926, 11124,
  23293, 51426, 832, 7068, 57646, 152002, 7531, 1398, 100101267, 221937,
  26873, 3709, 10576, 27040, 201176, 284403, 307, 6525, 5387, 1808,
  114907, 1952, 4091, 5096, 23451, 51207, 142891, 4084, 3797, 1185,
  5334, 217, 2042, 7781, 253461, 152404, 2992, 153478, 6695, 26249,
  662, 23036, 22852, 23411, 23387, 6786, 7690, 10021, 8925, 5595,
  63892, 23466, 11149, 139411, 755, 55347, 22820, 55111, 84444, 6258,
  57623, 1387, 3762, 5861, 3339, 114822, 129787, 196528, 54437, 10395,
  4023, 5095, 894, 3156, 1877, 283871, 2673, 23157, 10512, 2258,
  79750, 84665, 26091, 5978, 8751, 79695, 2768, 155382, 84163, 7709,
  4092, 1837, 2064, 4629, 55122, 8882, 150962, 23013, 8742, 489,
  201134, 55114, 84264, 64428, 4089, 861, 7422
) %>% as.character()

# ----------------- Load CorMotif Groups -----------------
# 0.1 µM
prob_1_0.1 <- as.character(read.csv("data/prob_1_0.1.csv")$Entrez_ID)
prob_2_0.1 <- as.character(read.csv("data/prob_2_0.1.csv")$Entrez_ID)
prob_3_0.1 <- as.character(read.csv("data/prob_3_0.1.csv")$Entrez_ID)

# 0.5 µM
prob_1_0.5 <- as.character(read.csv("data/prob_1_0.5.csv")$Entrez_ID)
prob_2_0.5 <- as.character(read.csv("data/prob_2_0.5.csv")$Entrez_ID)
prob_3_0.5 <- as.character(read.csv("data/prob_3_0.5.csv")$Entrez_ID)
prob_4_0.5 <- as.character(read.csv("data/prob_4_0.5.csv")$Entrez_ID)
prob_5_0.5 <- as.character(read.csv("data/prob_5_0.5.csv")$Entrez_ID)

# ----------------- Annotate Response Groups -----------------
df_0.1 <- data.frame(Entrez_ID = unique(c(prob_1_0.1, prob_2_0.1, prob_3_0.1))) %>%
  filter(Entrez_ID %in% entrez_ids) %>%
  mutate(
    Response_Group = case_when(
      Entrez_ID %in% prob_1_0.1 ~ "Non response",
      Entrez_ID %in% prob_2_0.1 ~ "CX-DOX mid-late response",
      Entrez_ID %in% prob_3_0.1 ~ "DOX only mid-late response"
    ),
    Concentration = "0.1 µM"
  )

df_0.5 <- data.frame(Entrez_ID = unique(c(prob_1_0.5, prob_2_0.5, prob_3_0.5, prob_4_0.5, prob_5_0.5))) %>%
  filter(Entrez_ID %in% entrez_ids) %>%
  mutate(
    Response_Group = case_when(
      Entrez_ID %in% prob_1_0.5 ~ "Non response",
      Entrez_ID %in% prob_2_0.5 ~ "DOX-specific response",
      Entrez_ID %in% prob_3_0.5 ~ "DOX only mid-late response",
      Entrez_ID %in% prob_4_0.5 ~ "CX total + DOX early response",
      Entrez_ID %in% prob_5_0.5 ~ "DOX early + CX-DOX mid-late response"
    ),
    Concentration = "0.5 µM"
  )

df_combined <- bind_rows(df_0.1, df_0.5)

# ----------------- Load and Aggregate DEG Data -----------------
deg_files <- list(
  "CX_0.1_3" = "data/DEGs/Toptable_CX_0.1_3.csv",
  "CX_0.1_24" = "data/DEGs/Toptable_CX_0.1_24.csv",
  "CX_0.1_48" = "data/DEGs/Toptable_CX_0.1_48.csv",
  "CX_0.5_3" = "data/DEGs/Toptable_CX_0.5_3.csv",
  "CX_0.5_24" = "data/DEGs/Toptable_CX_0.5_24.csv",
  "CX_0.5_48" = "data/DEGs/Toptable_CX_0.5_48.csv",
  "DOX_0.1_3" = "data/DEGs/Toptable_DOX_0.1_3.csv",
  "DOX_0.1_24" = "data/DEGs/Toptable_DOX_0.1_24.csv",
  "DOX_0.1_48" = "data/DEGs/Toptable_DOX_0.1_48.csv",
  "DOX_0.5_3" = "data/DEGs/Toptable_DOX_0.5_3.csv",
  "DOX_0.5_24" = "data/DEGs/Toptable_DOX_0.5_24.csv",
  "DOX_0.5_48" = "data/DEGs/Toptable_DOX_0.5_48.csv"
)

deg_all <- map2_dfr(deg_files, names(deg_files), ~{
  read_csv(.x, show_col_types = FALSE) %>% mutate(Condition = .y)
}) %>%
  mutate(
    Entrez_ID = as.character(Entrez_ID),
    Concentration = case_when(
      str_detect(Condition, "0.1") ~ "0.1 µM",
      str_detect(Condition, "0.5") ~ "0.5 µM"
    )
  )

deg_mean <- deg_all %>%
  group_by(Entrez_ID, Concentration) %>%
  summarise(logFC = mean(logFC, na.rm = TRUE), .groups = "drop")

# ----------------- Merge and Annotate -----------------
af_logfc <- inner_join(df_combined, deg_mean, by = c("Entrez_ID", "Concentration")) %>%
  mutate(
    Gene = mapIds(org.Hs.eg.db, keys = Entrez_ID,
                  column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"),
    Group = paste(Response_Group, Concentration)
  )

# ----------------- Response Group Order -----------------
group_levels <- c(
  "Non response",
  "CX-DOX mid-late response",
  "DOX-specific response",
  "DOX only mid-late response",
  "CX total + DOX early response",
  "DOX early + CX-DOX mid-late response"
)
af_logfc$Response_Group <- factor(af_logfc$Response_Group, levels = group_levels)

# ----------------- Wilcoxon Test -----------------
comparison_map <- list(
  "CX-DOX mid-late response 0.1 µM" = "Non response 0.1 µM",
  "DOX only mid-late response 0.1 µM" = "Non response 0.1 µM",
  "DOX-specific response 0.5 µM" = "Non response 0.5 µM",
  "DOX only mid-late response 0.5 µM" = "Non response 0.5 µM",
  "CX total + DOX early response 0.5 µM" = "Non response 0.5 µM",
  "DOX early + CX-DOX mid-late response 0.5 µM" = "Non response 0.5 µM"
)

wilcoxon_results <- lapply(names(comparison_map), function(resp_group) {
  control_group <- comparison_map[[resp_group]]
  resp_vals <- af_logfc$logFC[af_logfc$Group == resp_group]
  control_vals <- af_logfc$logFC[af_logfc$Group == control_group]

  if (length(resp_vals) >= 2 && length(control_vals) >= 2) {
    test_result <- wilcox.test(resp_vals, control_vals)
    pval <- test_result$p.value
    label <- case_when(
      pval < 0.001 ~ "***",
      pval < 0.01  ~ "**",
      pval < 0.05  ~ "*",
      TRUE ~ ""
    )
    y_pos <- max(resp_vals, na.rm = TRUE) + 0.5
    return(data.frame(Group = resp_group, y_pos = y_pos, label = label, P_Value = signif(pval, 4)))
  }
  return(NULL)
}) %>% bind_rows()

label_data <- wilcoxon_results %>%
  separate(Group, into = c("Response_Group", "Concentration"), sep = " (?=0\\.)") %>%
  mutate(Response_Group = factor(Response_Group, levels = group_levels))

# ----------------- Plot -----------------
ggplot(af_logfc, aes(x = Response_Group, y = logFC, fill = Response_Group)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7) +
  geom_text(data = label_data, aes(x = Response_Group, y = y_pos, label = label),
            inherit.aes = FALSE, size = 5, color = "black") +
  facet_wrap(~Concentration, scales = "free_x") +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Mean logFC of Atrial Fibrillation Genes\nAcross CorMotif Response Groups",
    x = "Response Group",
    y = "Mean Log Fold Change",
    fill = "Response Group"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 13),
    strip.text = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)
  )




############ HEART FAILURE GENE ANALYSIS #######################

# ----------------- Load Libraries -----------------
library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)

# ----------------- Heart Failure Entrez IDs -----------------
entrez_ids <- c(
  9709, 8882, 4023, 29959, 5496, 3992, 9415, 5308, 1026, 54437, 79068, 10221,
  9031, 1187, 1952, 3705, 84722, 7273, 23293, 155382, 9531, 602, 27258, 84163,
  81846, 79933, 56911, 64753, 93210, 1021, 283450, 5998, 57602, 114991, 7073,
  3156, 100101267, 22996, 285025, 11080, 11124, 54810, 7531, 27241, 4774, 57794,
  463, 91319, 6598, 9640, 2186, 26010, 80816, 571, 88, 51652, 64788, 90523, 2969,
  7781, 80777, 10725, 23387, 817, 134728, 8842, 949, 6934, 129787, 10327, 202052,
  2318, 5578, 6801, 6311, 10019, 80724, 217, 84909, 388591, 55101, 9839, 27161,
  5310, 387119, 4641, 5587, 55188, 222553, 9960, 22852, 10087, 9570, 54497,
  200942, 26249, 4137, 375056, 5409, 64116, 8291, 22876, 339855, 4864, 5142,
  221692, 55023, 51426, 6146, 84251, 8189, 27332, 57099, 1869, 1112, 23327,
  11264, 6001
) %>% as.character()

# ----------------- Load CorMotif Groups -----------------
# 0.1 µM
prob_1_0.1 <- as.character(read.csv("data/prob_1_0.1.csv")$Entrez_ID)
prob_2_0.1 <- as.character(read.csv("data/prob_2_0.1.csv")$Entrez_ID)
prob_3_0.1 <- as.character(read.csv("data/prob_3_0.1.csv")$Entrez_ID)

# 0.5 µM
prob_1_0.5 <- as.character(read.csv("data/prob_1_0.5.csv")$Entrez_ID)
prob_2_0.5 <- as.character(read.csv("data/prob_2_0.5.csv")$Entrez_ID)
prob_3_0.5 <- as.character(read.csv("data/prob_3_0.5.csv")$Entrez_ID)
prob_4_0.5 <- as.character(read.csv("data/prob_4_0.5.csv")$Entrez_ID)
prob_5_0.5 <- as.character(read.csv("data/prob_5_0.5.csv")$Entrez_ID)

# ----------------- Annotate Response Groups -----------------
df_0.1 <- data.frame(Entrez_ID = unique(c(prob_1_0.1, prob_2_0.1, prob_3_0.1))) %>%
  filter(Entrez_ID %in% entrez_ids) %>%
  mutate(
    Response_Group = case_when(
      Entrez_ID %in% prob_1_0.1 ~ "Non response",
      Entrez_ID %in% prob_2_0.1 ~ "CX-DOX mid-late response",
      Entrez_ID %in% prob_3_0.1 ~ "DOX only mid-late response"
    ),
    Concentration = "0.1 µM"
  )

df_0.5 <- data.frame(Entrez_ID = unique(c(prob_1_0.5, prob_2_0.5, prob_3_0.5, prob_4_0.5, prob_5_0.5))) %>%
  filter(Entrez_ID %in% entrez_ids) %>%
  mutate(
    Response_Group = case_when(
      Entrez_ID %in% prob_1_0.5 ~ "Non response",
      Entrez_ID %in% prob_2_0.5 ~ "DOX-specific response",
      Entrez_ID %in% prob_3_0.5 ~ "DOX only mid-late response",
      Entrez_ID %in% prob_4_0.5 ~ "CX total + DOX early response",
      Entrez_ID %in% prob_5_0.5 ~ "DOX early + CX-DOX mid-late response"
    ),
    Concentration = "0.5 µM"
  )

df_combined <- bind_rows(df_0.1, df_0.5)

# ----------------- Load and Aggregate DEG Data -----------------
deg_files <- list(
  "CX_0.1_3" = "data/DEGs/Toptable_CX_0.1_3.csv",
  "CX_0.1_24" = "data/DEGs/Toptable_CX_0.1_24.csv",
  "CX_0.1_48" = "data/DEGs/Toptable_CX_0.1_48.csv",
  "CX_0.5_3" = "data/DEGs/Toptable_CX_0.5_3.csv",
  "CX_0.5_24" = "data/DEGs/Toptable_CX_0.5_24.csv",
  "CX_0.5_48" = "data/DEGs/Toptable_CX_0.5_48.csv",
  "DOX_0.1_3" = "data/DEGs/Toptable_DOX_0.1_3.csv",
  "DOX_0.1_24" = "data/DEGs/Toptable_DOX_0.1_24.csv",
  "DOX_0.1_48" = "data/DEGs/Toptable_DOX_0.1_48.csv",
  "DOX_0.5_3" = "data/DEGs/Toptable_DOX_0.5_3.csv",
  "DOX_0.5_24" = "data/DEGs/Toptable_DOX_0.5_24.csv",
  "DOX_0.5_48" = "data/DEGs/Toptable_DOX_0.5_48.csv"
)

deg_all <- map2_dfr(deg_files, names(deg_files), ~{
  read_csv(.x, show_col_types = FALSE) %>%
    mutate(Condition = .y)
}) %>%
  mutate(
    Entrez_ID = as.character(Entrez_ID),
    Concentration = case_when(
      str_detect(Condition, "0.1") ~ "0.1 µM",
      str_detect(Condition, "0.5") ~ "0.5 µM"
    )
  )

deg_mean <- deg_all %>%
  group_by(Entrez_ID, Concentration) %>%
  summarise(logFC = mean(logFC, na.rm = TRUE), .groups = "drop")

# ----------------- Merge and Annotate -----------------
hf_logfc <- inner_join(df_combined, deg_mean, by = c("Entrez_ID", "Concentration")) %>%
  mutate(
    Gene = mapIds(org.Hs.eg.db, keys = Entrez_ID,
                  column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"),
    Group = paste(Response_Group, Concentration)
  )

# ----------------- Define Group Order -----------------
group_levels <- c(
  "Non response",
  "CX-DOX mid-late response",
  "DOX-specific response",
  "DOX only mid-late response",
  "CX total + DOX early response",
  "DOX early + CX-DOX mid-late response"
)
hf_logfc$Response_Group <- factor(hf_logfc$Response_Group, levels = group_levels)

# ----------------- Wilcoxon Test -----------------
comparison_map <- list(
  "CX-DOX mid-late response 0.1 µM" = "Non response 0.1 µM",
  "DOX only mid-late response 0.1 µM" = "Non response 0.1 µM",
  "DOX-specific response 0.5 µM" = "Non response 0.5 µM",
  "DOX only mid-late response 0.5 µM" = "Non response 0.5 µM",
  "CX total + DOX early response 0.5 µM" = "Non response 0.5 µM",
  "DOX early + CX-DOX mid-late response 0.5 µM" = "Non response 0.5 µM"
)

wilcoxon_results <- lapply(names(comparison_map), function(resp_group) {
  control_group <- comparison_map[[resp_group]]
  resp_vals <- hf_logfc$logFC[hf_logfc$Group == resp_group]
  control_vals <- hf_logfc$logFC[hf_logfc$Group == control_group]

  if (length(resp_vals) >= 2 && length(control_vals) >= 2) {
    test_result <- wilcox.test(resp_vals, control_vals)
    pval <- test_result$p.value
    label <- case_when(
      pval < 0.001 ~ "***",
      pval < 0.01  ~ "**",
      pval < 0.05  ~ "*",
      TRUE ~ ""
    )
    y_pos <- max(resp_vals, na.rm = TRUE) + 0.5
    return(data.frame(Group = resp_group, y_pos = y_pos, label = label, P_Value = signif(pval, 4)))
  }
  return(NULL)
}) %>% bind_rows()

label_data <- wilcoxon_results %>%
  separate(Group, into = c("Response_Group", "Concentration"), sep = " (?=0\\.)") %>%
  mutate(Response_Group = factor(Response_Group, levels = group_levels))

# ----------------- Plot -----------------
ggplot(hf_logfc, aes(x = Response_Group, y = logFC, fill = Response_Group)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7) +
  geom_text(data = label_data, aes(x = Response_Group, y = y_pos, label = label),
            inherit.aes = FALSE, size = 5, color = "black") +
  facet_wrap(~Concentration, scales = "free_x") +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Mean logFC of Heart Failure Genes\nAcross CorMotif Response Groups",
    x = "Response Group",
    y = "Mean Log Fold Change",
    fill = "Response Group"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 13),
    strip.text = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)
  )
