library(DESeq2)
library(stringr)
library(magick)
library(rsvg)
library(RColorBrewer)
# Exp data
load("~/Documents/MUNKA/PhD/Single_cell/sequencing_results/Final_analysis_2022/data/Combined_DESeqData.RData")
TissueSpecFBDevelopment <- unique(df_combined$sample[df_combined$seqtype %in% "TissuePairedEnd"])
SimpleSugarSolid <- c("Xyl_S", "Ara_S", "Fru_S", "Man_S", "GA_S", "SoA_S", "Rib_S", "Con_S", "Glu_S")
SimpleSugarLiq <- c( "Gluco_P", "Copciglucose_S", "Fucose_P", "Galose_P", "GA_P", "Glyc_P", "Lac_P",
                     "Malt_P", "Mannit_P", "Rham_P", "Sorb_P", "Treh_P")
PolysaccharideSolid <- c("PGA_S","Rha_S", "Pec_S")
PolysaccharideLiq <- c("PGA_P", "Mannan_P", "Inulin_P", "GG_P", "AP_P", "Amy_P", "Ara_P",
                       "Gala_P", "GM_P", "Glucan_P", "Xylan_P", "XG_P", "AG_P", "Pectin_P", "RhaG_P",
                       "CopciCellobiose_S", "Copcicellulose_S", "CopciLignin_S",  "CopciPectine_S",  "CopciXylose_S")
LignocelluloseLiq <- c("Lime_P", "Poplar_P", "Typha_P", "Salix_P", "Reed_P",
                       "TM_P", "Ecane_P", "Alga_P", "Myc_P", "VA_P", "FA_P",
                       "CopciWheatbran_S", "CopciHorsemanure_S", "CopciOakleaves_S",
                       "CopciHay_S", "CopciCornstalk_S", "CopciApplepeels_S")
NitrogenLiq <- c("NHNO_P", "Pro_P",  "Trp_P", "Ile_P", "Arg_P", "Meti_P", "NaNO_P")
Others <- c("BSA_P", "BSA_S", "NoC_S", "NoN_S", "NoP_S", "NoS_S" )
df_combined$groups[df_combined$sample %in% SimpleSugarSolid] <- "SimpleSugarSolid"
df_combined$groups[df_combined$sample %in% SimpleSugarLiq] <- "SimpleSugarLiq"
df_combined$groups[df_combined$sample %in% PolysaccharideSolid] <- "PolysaccharideSolid"
df_combined$groups[df_combined$sample %in% PolysaccharideLiq] <- "PolysaccharideLiq"
df_combined$groups[df_combined$sample %in% LignocelluloseLiq] <- "LignocelluloseLiq"
df_combined$groups[df_combined$sample %in% NitrogenLiq] <- "NitrogenLiq"
df_combined$groups[df_combined$sample %in% Others] <- "Others"
df_combined$groups[df_combined$sample %in% TissueSpecFBDevelopment] <- "TissueSpecFBDevelopment"
df_combined$sample <- str_remove(df_combined$sample, "Copci")
SimpleSugarLiq <-  str_remove(SimpleSugarLiq, "Copci")
PolysaccharideLiq <-  str_remove(PolysaccharideLiq, "Copci")
LignocelluloseLiq <-  str_remove(LignocelluloseLiq, "Copci")
levels_order <- c("VM_VM", "H1_VM", "H1_cap","H2_Nodulus", "H3_Nodulus",  "P0_Nodulus", "P1_Nodulus",  "P2_Nodulus",
                  "P0_External_nodulus", "P1_External_nodulus", "P2_External_nodulus",
                  "P1_Stipe", "P2_Stipe", "P2_Central_stipe",
                  "H2_UV", "H3_UV", "P0_UV",  "P1_UV", "P2_UV",
                  "H3_Cap",  "P0_Cap",   "P1_Cap", "P2_Cap", "P1_Upper_cap", "P2_Upper_cap",
                  "P1_PV", "P2_PV",  "P1_Gill",  "P2_Gill",
                  SimpleSugarSolid, SimpleSugarLiq, PolysaccharideSolid, PolysaccharideLiq, LignocelluloseLiq,
                  NitrogenLiq, Others)
levels_order[duplicated(levels_order)]
df_combined$order <- factor(df_combined$sample, levels =  levels_order)
df_combined$plot_facets[df_combined$sample %in% c(SimpleSugarSolid, SimpleSugarLiq)] <- "Simple_sugars"
df_combined$plot_facets[df_combined$sample %in% c(PolysaccharideSolid, PolysaccharideLiq)] <- "Polysaccharide"
df_combined$plot_facets[df_combined$sample %in% c(LignocelluloseLiq, NitrogenLiq, Others)] <- "Lignocellulose_Nitrogen_etc"
df_combined$plot_facets[df_combined$sample %in% c(TissueSpecFBDevelopment)] <- "TissueSpecFBDevelopment"

df_combined$plot_facets <- factor(df_combined$plot_facets, levels = c("TissueSpecFBDevelopment",
                                                                      "Simple_sugars", "Polysaccharide",
                                                                      "Lignocellulose_Nitrogen_etc"))

df_combined$fill <- str_remove(df_combined$sample, "([^_]+)_")
df_combined$fill[df_combined$sample %in% "H1_cap"] <- "H1"
df_combined$fill[df_combined$sample %in% "H1_VM"] <- "H1"
df_combined$fill[!df_combined$groups %in% "TissueSpecFBDevelopment"] <- "VM"
df_combined$fill <- factor(df_combined$fill, levels = c("VM", "H1", "Nodulus",
                                                        "External_nodulus",  "Stipe", "Central_stipe",
                                                        "UV", "Cap", "Upper_cap", "PV", "Gill"))
ddSE <- DESeqDataSetFromMatrix(countData = counts_combined, colData = df_combined, design = ~ sample)
normalizationFactors(ddSE) <- EDASeqNormFactors_Combined
vsd_smoc2 <- vst(ddSE, blind = TRUE)
exp_data <- assay(vsd_smoc2)
usethis::use_data(exp_data, df_combined, internal = TRUE)
# Image data
load("./R/sysdata.rda")
df_plot <- df_combined[df_combined$groups %in% "TissueSpecFBDevelopment", ]
samples <- unique(df_plot$sample)
df_plot <- data.frame(images = sapply(str_split(samples, "_"), function(x) x[1]),
                      areas = samples, coords = rep(NA, length(samples)),
                      stringsAsFactors = FALSE)
df_plot <- df_plot[order(df_plot$images, df_plot$areas),]
df_plot <- df_plot[c(29, 1:28),]
### Define tissue_names, that have multiple areas
df_plot$areas[str_detect(df_plot$areas, "External")] <- str_c(df_plot$areas[str_detect(df_plot$areas, "External")], "_A")
df_temp <- df_plot[str_detect(df_plot$areas, "External"),]
df_temp$areas <- str_replace(df_temp$areas, "_A", "_B")
df_plot <- rbind(df_plot, df_temp)
#### P0 has only 1 external nodulus sample
df_plot <- df_plot[-30,]
df_plot$areas[df_plot$areas %in% "P0_External_nodulus_A"] <- "P0_External_nodulus"
df_plot$areas[str_detect(df_plot$areas, "Gill")] <- str_c(df_plot$areas[str_detect(df_plot$areas, "Gill")], "_A")
df_temp <- df_plot[str_detect(df_plot$areas, "Gill"),]
df_temp$areas <- str_replace(df_temp$areas, "_A", "_B")
df_plot <- rbind(df_plot, df_temp)
df_plot$areas[str_detect(df_plot$areas, "PV")] <- str_c(df_plot$areas[str_detect(df_plot$areas, "PV")], "_A")
df_temp <- df_plot[str_detect(df_plot$areas, "PV"),]
df_temp$areas <- str_replace(df_temp$areas, "_A", "_B")
df_plot <- rbind(df_plot, df_temp)
df_plot <- df_plot[order(df_plot$images, df_plot$areas),]
df_plot <- df_plot[c(35, 1:34),]
df_plot <- df_plot[,1:2] # only the first two columns are needed

image_data <- image_data_process(df_plot)
image_df <- image_data$image_df


### mean expression data
df_combined2 <- df_combined[df_combined$groups %in% "TissueSpecFBDevelopment", ]
expression_data <- exp_data[,colnames(exp_data) %in% df_combined2$Library_ID_final]
mean_exp_data <- mean_expression(expression_data = expression_data,
                                 coldata = df_combined2,
                                 sample_col = "sample",
                                 id_col = "Library_ID_final")

usethis::use_data(exp_data, df_combined, image_df, mean_exp_data, internal = TRUE, overwrite = TRUE)

# Test copci_image_data and copci plot

image_data <- copci_image_data()
copci_dev_heatmap(gene = "CopciAB_411205", image_data = image_data)
