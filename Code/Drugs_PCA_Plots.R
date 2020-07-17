suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(pcaMethods))

# PCA - Positive Mode Data ------------------------------------------------

# input data w/ metadata appended
df <- fread("input/20200225_DarkWebDrug_Pos.csv", sep=",", header=TRUE)

# subset data to only include samples
df <- subset(df, df$ATTRIBUTE_SampleType == "reference material")

# calculate
mypca <- pca(df[,-c(1:53)], method="nipals", center=TRUE, nPcs = 9, scale = "pareto")
df <- cbind(scores(mypca),df)
  positive_scores <- df
loadings <- as.data.frame(loadings(mypca))
rownames <- rownames(loadings)
loadings <- cbind(rownames,loadings)
rownames(loadings) <- NULL
test_loading_names <- loadings
test_loading_names <- cbind(substr(test_loading_names$rownames, 0, 9),loadings)
colnames(test_loading_names)[1] <- "mz"
test_loading_names <- separate(test_loading_names, col="mz", remove= TRUE, "_", into=c("feature_number","mz"))

# add annotations
GNPS_hits <- fread("input/GNPS_libraryhits/FEATURE-BASED-MOLECULAR-NETWORKING-0834479d-view_all_annotations_DB-main.tsv")
test_loading_names <- merge(test_loading_names, GNPS_hits[,c("#Scan#","Compound_Name")], by.x="feature_number", by.y="#Scan#", all = TRUE)

# ordering manufacturer info
df$ATTRIBUTE_Info_Label_Manufacturer <- factor(df$ATTRIBUTE_Info_Label_Manufacturer, 
                                               levels = c('CIPLA LTD','Combitic Global Caplet PVT. LTD.','Healing Pharma',
                                                          'fortune health care','Laborate Pharmaceuticals India LTD.','pharose remedies Ltd.',
                                                          'Naman Pharma Drugs','Signature Phytochemical Industries',
                                                          'Skohind Labs strategic alliance with Rousel Laboratories',
                                                          'uni medicolabs','Medreich limited','Centurion Laboratories PVT. LTD.'))
index_man <- fread("input/renumbering_manuf_index.txt", sep="\t", header=TRUE)
df <- merge(index_man, df, by.x="manufacturer", by.y="ATTRIBUTE_Info_Label_Manufacturer")

# Plots
plot_value <- df$ATTRIBUTE_API

score_plot <- ggplot(df, aes(as.numeric(PC2), as.numeric(PC3)))+
  geom_point(aes(colour=as.factor(plot_value)), pch=16, cex=1.0, alpha=0.7) +
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  #geom_text(aes(label=filename), size = 2)+
  theme_minimal()+
  theme(panel.background = element_rect(fill="grey85"),
        panel.grid.major=element_line(colour="grey90",size=0.25, linetype="dashed"),
        panel.grid.minor=element_line(colour="grey90",size=0.25, linetype="dashed"),
        panel.border = element_rect(colour ="black", size=0.25, linetype = "solid", fill = NA),
        axis.ticks=element_line(colour ="black",size=0.25, linetype="solid"),
        axis.text=element_text(colour="black",size=6),
        axis.text.x=element_text(colour="black",size=6),
        axis.line=element_line(colour="black",size=0.25, linetype="solid"),                  
        axis.title=element_text(colour="black",size=6),
        strip.text.x = element_text(colour="black",size=6),
        aspect.ratio=1,
        legend.position="none",
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size=6)) +
  labs(x=paste("PC2 (",round(mypca@R2[2] *100,2), "%)"), y=paste("PC3 (",round(mypca@R2[3] *100,2), "%)"))
print(score_plot)
ggsave(score_plot, filename="R_output/Drugs_Pos_PCA_PC2vsPC3_Score_All_API.pdf", width = 2, height = 2, units = "in", dpi=300, useDingbats=FALSE)

# Plots
plot_value <- df$ATTRIBUTE_API_Drug2

score_plot <- ggplot(df, aes(as.numeric(PC2), as.numeric(PC3)))+
  geom_point(aes(colour=as.factor(plot_value)), pch=16, cex=1.0, alpha=0.7) +
  scale_color_brewer(palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  #geom_text(aes(label=filename), size = 2)+
  theme_minimal()+
  theme(panel.background = element_rect(fill="grey85"),
        panel.grid.major=element_line(colour="grey90",size=0.25, linetype="dashed"),
        panel.grid.minor=element_line(colour="grey90",size=0.25, linetype="dashed"),
        panel.border = element_rect(colour ="black", size=0.25, linetype = "solid", fill = NA),
        axis.ticks=element_line(colour ="black",size=0.25, linetype="solid"),
        axis.text=element_text(colour="black",size=6),
        axis.text.x=element_text(colour="black",size=6),
        axis.line=element_line(colour="black",size=0.25, linetype="solid"),                  
        axis.title=element_text(colour="black",size=6),
        strip.text.x = element_text(colour="black",size=6),
        aspect.ratio=1,
        legend.position="none",
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size=6)) +
  labs(x=paste("PC2 (",round(mypca@R2[2] *100,2), "%)"), y=paste("PC3 (",round(mypca@R2[3] *100,2), "%)"))
print(score_plot)
ggsave(score_plot, filename="R_output/Drugs_Pos_PCA_PC2vsPC3_Score_All_APIDrug2.pdf", width = 2, height = 2, units = "in", dpi=300, useDingbats=FALSE)

loading_plot <- ggplot(test_loading_names, aes(PC2, PC3, label=feature_number))+
  geom_point(colour="#525252", pch=16,cex=0.5, alpha=0.6)+
  geom_text(check_overlap = TRUE, size = 2, vjust = -1) +
  geom_point(data= subset(test_loading_names, test_loading_names$Compound_Name != "NA"), 
             aes(PC2, PC3), colour="red", pch=16,cex=1.0, alpha=0.6)+
  theme_minimal()+
  theme(panel.grid.major=element_line(colour="grey90",size=0.25, linetype="dashed"),
        panel.grid.minor=element_line(colour="grey90",size=0.25, linetype="dashed"),
        panel.border = element_rect(colour ="black", size=0.25, linetype = "solid", fill = NA),
        axis.ticks=element_line(colour ="black",size=0.25, linetype="solid"),
        axis.text=element_text(colour="black",size=6),
        axis.text.x=element_text(colour="black",size=6),
        axis.line=element_line(colour="black",size=0.25, linetype="solid"),                  
        axis.title=element_text(colour="black",size=6),
        strip.text.x = element_text(colour="black",size=6),
        aspect.ratio=1,
        legend.position="none",
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size=6)) +
  labs(x="PC2 loadings", y="PC3 loadings")
print(loading_plot)
ggsave(loading_plot, filename="R_output/Drugs_Pos_PCA_PC2vsPC3_Loading.pdf", width = 3.25, height = 3.25, units = "in", dpi=300, useDingbats=FALSE)

# facet on Manufacturer and Package ID
score_plot <- ggplot(df, aes(as.numeric(PC2), as.numeric(PC3)))+
  geom_point(data = select(df, -c(number_man,Cbind_Sample_ID)), aes(colour=ATTRIBUTE_API), pch=16, cex=0.5, alpha=0.5) +
  geom_point(aes(fill=ATTRIBUTE_API), colour="black", pch=21, cex=1.0, alpha=1.0) +
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  theme_minimal()+
  facet_wrap(~number_man~Cbind_Sample_ID, ncol = 7)+
  theme(panel.background = element_rect(fill="grey85"),
        panel.grid.major=element_line(colour="grey90",size=0.25, linetype="dashed"),
        panel.grid.minor=element_line(colour="grey90",size=0.25, linetype="dashed"),
        panel.border = element_rect(colour ="black", size=0.25, linetype = "solid", fill = NA),
        axis.ticks=element_line(colour ="black",size=0.25, linetype="solid"),
        axis.text=element_text(colour="black",size=6),
        axis.text.x=element_text(colour="black",size=6),
        axis.line=element_line(colour="black",size=0.25, linetype="solid"),                  
        axis.title=element_text(colour="black",size=6),
        strip.text.x = element_text(colour="black",size=6),
        aspect.ratio=1,
        legend.position="none",
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size=6)) +
  labs(x=paste("PC2 (",round(mypca@R2[2] *100,2), "%)"), y=paste("PC3 (",round(mypca@R2[3] *100,2), "%)"))
print(score_plot)
ggsave(score_plot, filename="R_output/Drugs_Pos_PCA_PC2vsPC3_Score_All_API_Facet.pdf", width = 7.1, height = 11, units = "in", dpi=300, useDingbats=FALSE)

# variance plots
variance <- round(mypca@R2*100,2)
cumul_variance <- round(mypca@R2cum*100,2)
PCs_num <- c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9")
plot_variance <- data.frame(PCs_num, variance, cumul_variance)

# facet on Manufacturer and Package ID
variance <- ggplot(plot_variance, aes(as.numeric(PCs_num), as.numeric(variance)))+
  geom_point(colour="black", pch=16, cex=1.25, alpha=1.0)+
  geom_path(aes(as.numeric(PCs_num), as.numeric(variance)),colour="black")+
  geom_point(aes(x=as.numeric(PCs_num), y=as.numeric(cumul_variance)), colour="red", pch=1, cex=1.25, alpha=1.0)+
  geom_path(aes(as.numeric(PCs_num), as.numeric(cumul_variance)),colour="red")+
  theme_minimal()+
  theme(panel.background = element_rect(fill="grey85"),
        panel.grid.major=element_line(colour="grey90",size=0.25, linetype="dashed"),
        panel.grid.minor=element_line(colour="grey90",size=0.25, linetype="dashed"),
        panel.border = element_rect(colour ="black", size=0.25, linetype = "solid", fill = NA),
        axis.ticks=element_line(colour ="black",size=0.25, linetype="solid"),
        axis.text=element_text(colour="black",size=6),
        axis.text.x=element_text(colour="black",size=6),
        axis.line=element_line(colour="black",size=0.25, linetype="solid"),                  
        axis.title=element_text(colour="black",size=6),
        strip.text.x = element_text(colour="black",size=6),
        aspect.ratio=1,
        legend.position="none",
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size=6)) +
  labs(x="PCs", y="% variance")
print(variance)
ggsave(variance, filename="R_output/Drugs_Pos_PCA_Variance.pdf", width = 2.25, height = 2.25, units = "in", dpi=300, useDingbats=FALSE)


# PCA - Negative Mode Data ------------------------------------------------

# input data w/ metadata appended
df <- fread("input/20200225_DarkWebDrug_Neg.csv", sep=",", header=TRUE)

# subset data to only include samples
df <- subset(df, df$ATTRIBUTE_SampleType == "reference material")

# calculate
mypca <- pca(df[,-c(1:53)], method="nipals", center=TRUE, nPcs = 9, scale = "pareto")
df <- cbind(scores(mypca),df)
  negative_scores <- df
loadings <- as.data.frame(loadings(mypca))
rownames <- rownames(loadings)
loadings <- cbind(rownames,loadings)
rownames(loadings) <- NULL
test_loading_names <- loadings
test_loading_names <- cbind(substr(test_loading_names$rownames, 0, 9),loadings)
colnames(test_loading_names)[1] <- "mz"
test_loading_names <- separate(test_loading_names, col="mz", remove= TRUE, "_", into=c("feature_number","mz"))

# add annotations
GNPS_hits <- fread("input/GNPS_libraryhits/FEATURE-BASED-MOLECULAR-NETWORKING-a40b5a63-view_all_annotations_DB-main.tsv")
test_loading_names <- merge(test_loading_names, GNPS_hits[,c("#Scan#","Compound_Name")], by.x="feature_number", by.y="#Scan#", all = TRUE)

# ordering manufacturer info
df$ATTRIBUTE_Info_Label_Manufacturer <- factor(df$ATTRIBUTE_Info_Label_Manufacturer, 
                                               levels = c('CIPLA LTD','Combitic Global Caplet PVT. LTD.','Healing Pharma',
                                                          'fortune health care','Laborate Pharmaceuticals India LTD.','pharose remedies Ltd.',
                                                          'Naman Pharma Drugs','Signature Phytochemical Industries',
                                                          'Skohind Labs strategic alliance with Rousel Laboratories',
                                                          'uni medicolabs','Medreich limited','Centurion Laboratories PVT. LTD.'))
index_man <- fread("input/renumbering_manuf_index.txt", sep="\t", header=TRUE)
df <- merge(index_man, df, by.x="manufacturer", by.y="ATTRIBUTE_Info_Label_Manufacturer")

# Plots
plot_value <- df$ATTRIBUTE_API

score_plot <- ggplot(df, aes(as.numeric(PC2), as.numeric(PC3)))+
  geom_point(aes(colour=as.factor(plot_value)), pch=16, cex=1.0, alpha=0.7) +
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  #geom_text(aes(label=filename), size = 2)+
  theme_minimal()+
  theme(panel.background = element_rect(fill="grey85"),
        panel.grid.major=element_line(colour="grey90",size=0.25, linetype="dashed"),
        panel.grid.minor=element_line(colour="grey90",size=0.25, linetype="dashed"),
        panel.border = element_rect(colour ="black", size=0.25, linetype = "solid", fill = NA),
        axis.ticks=element_line(colour ="black",size=0.25, linetype="solid"),
        axis.text=element_text(colour="black",size=6),
        axis.text.x=element_text(colour="black",size=6),
        axis.line=element_line(colour="black",size=0.25, linetype="solid"),                  
        axis.title=element_text(colour="black",size=6),
        strip.text.x = element_text(colour="black",size=6),
        aspect.ratio=1,
        legend.position="none",
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size=6)) +
  labs(x=paste("PC2 (",round(mypca@R2[2] *100,2), "%)"), y=paste("PC3 (",round(mypca@R2[3] *100,2), "%)"))
print(score_plot)
ggsave(score_plot, filename="R_output/Drugs_Neg_PCA_PC2vsPC3_Score_All_API.pdf", width = 2, height = 2, units = "in", dpi=300, useDingbats=FALSE)

# Plots
plot_value <- df$ATTRIBUTE_API_Drug2

score_plot <- ggplot(df, aes(as.numeric(PC2), as.numeric(PC3)))+
  geom_point(aes(colour=as.factor(plot_value)), pch=16, cex=1.0, alpha=0.7) +
  scale_color_brewer(palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  #geom_text(aes(label=filename), size = 2)+
  theme_minimal()+
  theme(panel.background = element_rect(fill="grey85"),
        panel.grid.major=element_line(colour="grey90",size=0.25, linetype="dashed"),
        panel.grid.minor=element_line(colour="grey90",size=0.25, linetype="dashed"),
        panel.border = element_rect(colour ="black", size=0.25, linetype = "solid", fill = NA),
        axis.ticks=element_line(colour ="black",size=0.25, linetype="solid"),
        axis.text=element_text(colour="black",size=6),
        axis.text.x=element_text(colour="black",size=6),
        axis.line=element_line(colour="black",size=0.25, linetype="solid"),                  
        axis.title=element_text(colour="black",size=6),
        strip.text.x = element_text(colour="black",size=6),
        aspect.ratio=1,
        legend.position="none",
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size=6)) +
  labs(x=paste("PC2 (",round(mypca@R2[2] *100,2), "%)"), y=paste("PC3 (",round(mypca@R2[3] *100,2), "%)"))
print(score_plot)
ggsave(score_plot, filename="R_output/Drugs_Neg_PCA_PC2vsPC3_Score_All_APIDrug2.pdf", width = 2, height = 2, units = "in", dpi=300, useDingbats=FALSE)

loading_plot <- ggplot(test_loading_names, aes(PC2, PC3, label=feature_number))+
  geom_point(colour="#525252", pch=16,cex=0.5, alpha=0.6)+
  geom_text(check_overlap = TRUE, size = 2, vjust = -1) +
  geom_point(data= subset(test_loading_names, test_loading_names$Compound_Name != "NA"), 
             aes(PC2, PC3), colour="red", pch=16,cex=1.0, alpha=0.6)+
  theme_minimal()+
  theme(panel.grid.major=element_line(colour="grey90",size=0.25, linetype="dashed"),
        panel.grid.minor=element_line(colour="grey90",size=0.25, linetype="dashed"),
        panel.border = element_rect(colour ="black", size=0.25, linetype = "solid", fill = NA),
        axis.ticks=element_line(colour ="black",size=0.25, linetype="solid"),
        axis.text=element_text(colour="black",size=6),
        axis.text.x=element_text(colour="black",size=6),
        axis.line=element_line(colour="black",size=0.25, linetype="solid"),                  
        axis.title=element_text(colour="black",size=6),
        strip.text.x = element_text(colour="black",size=6),
        aspect.ratio=1,
        legend.position="none",
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size=6)) +
  labs(x="PC2 loadings", y="PC3 loadings")
print(loading_plot)
ggsave(loading_plot, filename="R_output/Drugs_Neg_PCA_PC2vsPC3_Loading.pdf", width = 3.25, height = 3.25, units = "in", dpi=300, useDingbats=FALSE)

# facet on Manufacturer and Package ID
score_plot <- ggplot(df, aes(as.numeric(PC2), as.numeric(PC3)))+
  geom_point(data = select(df, -c(number_man,Cbind_Sample_ID)), aes(colour=ATTRIBUTE_API), pch=16, cex=0.5, alpha=0.5) +
  geom_point(aes(fill=ATTRIBUTE_API), colour="black", pch=21, cex=1.0, alpha=1.0) +
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  theme_minimal()+
  facet_wrap(~number_man~Cbind_Sample_ID, ncol = 7)+
  theme(panel.background = element_rect(fill="grey85"),
        panel.grid.major=element_line(colour="grey90",size=0.25, linetype="dashed"),
        panel.grid.minor=element_line(colour="grey90",size=0.25, linetype="dashed"),
        panel.border = element_rect(colour ="black", size=0.25, linetype = "solid", fill = NA),
        axis.ticks=element_line(colour ="black",size=0.25, linetype="solid"),
        axis.text=element_text(colour="black",size=6),
        axis.text.x=element_text(colour="black",size=6),
        axis.line=element_line(colour="black",size=0.25, linetype="solid"),                  
        axis.title=element_text(colour="black",size=6),
        strip.text.x = element_text(colour="black",size=6),
        aspect.ratio=1,
        legend.position="none",
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size=6)) +
  labs(x=paste("PC2 (",round(mypca@R2[2] *100,2), "%)"), y=paste("PC3 (",round(mypca@R2[3] *100,2), "%)"))
print(score_plot)
ggsave(score_plot, filename="R_output/Drugs_Neg_PCA_PC2vsPC3_Score_All_API_Facet.pdf", width = 7.1, height = 11, units = "in", dpi=300, useDingbats=FALSE)


# variance plots
variance <- round(mypca@R2*100,2)
cumul_variance <- round(mypca@R2cum*100,2)
PCs_num <- c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9")
plot_variance <- data.frame(PCs_num, variance, cumul_variance)

# facet on Manufacturer and Package ID
variance <- ggplot(plot_variance, aes(as.numeric(PCs_num), as.numeric(variance)))+
  geom_point(colour="black", pch=16, cex=1.25, alpha=1.0)+
  geom_path(aes(as.numeric(PCs_num), as.numeric(variance)),colour="black")+
  geom_point(aes(x=as.numeric(PCs_num), y=as.numeric(cumul_variance)), colour="red", pch=1, cex=1.25, alpha=1.0)+
  geom_path(aes(as.numeric(PCs_num), as.numeric(cumul_variance)),colour="red")+
  theme_minimal()+
  theme(panel.background = element_rect(fill="grey85"),
        panel.grid.major=element_line(colour="grey90",size=0.25, linetype="dashed"),
        panel.grid.minor=element_line(colour="grey90",size=0.25, linetype="dashed"),
        panel.border = element_rect(colour ="black", size=0.25, linetype = "solid", fill = NA),
        axis.ticks=element_line(colour ="black",size=0.25, linetype="solid"),
        axis.text=element_text(colour="black",size=6),
        axis.text.x=element_text(colour="black",size=6),
        axis.line=element_line(colour="black",size=0.25, linetype="solid"),                  
        axis.title=element_text(colour="black",size=6),
        strip.text.x = element_text(colour="black",size=6),
        aspect.ratio=1,
        legend.position="none",
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size=6)) +
  labs(x="PCs", y="% variance")
print(variance)
ggsave(variance, filename="R_output/Drugs_Neg_PCA_Variance.pdf", width = 2.25, height = 2.25, units = "in", dpi=300, useDingbats=FALSE)


# FUSED

negative_scores <- negative_scores %>% select(Cbind_Sample_ID,ATTRIBUTE_Rep_Day,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9)
positive_scores <- positive_scores %>% select(Cbind_Sample_ID,ATTRIBUTE_Rep_Day,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9)

fused_scores <- left_join(positive_scores, negative_scores, by=c("Cbind_Sample_ID","ATTRIBUTE_Rep_Day"))
fused_scores <- fused_scores[complete.cases(fused_scores), ]

# calculate
mypca <- pca(fused_scores[,-c(1:2)], method="nipals", center=TRUE, nPcs = 8, scale = "none")
df <- cbind(scores(mypca),fused_scores[,c(1:2)])

# get metadata
metadata <- fread("input/20200225_DarkWebDrug_Neg.csv", sep=",", header=TRUE)
metadata <- metadata %>% select(Cbind_Sample_ID,ATTRIBUTE_Rep_Day,ATTRIBUTE_Info_Label_Manufacturer,ATTRIBUTE_API,ATTRIBUTE_API_Drug2)


df <- left_join(metadata, df, by=c("Cbind_Sample_ID","ATTRIBUTE_Rep_Day"))
df <- subset(df, df$Cbind_Sample_ID != "not collected")


# ordering manufacturer info
df$ATTRIBUTE_Info_Label_Manufacturer <- factor(df$ATTRIBUTE_Info_Label_Manufacturer, 
                                               levels = c('CIPLA LTD','Combitic Global Caplet PVT. LTD.','Healing Pharma',
                                                          'fortune health care','Laborate Pharmaceuticals India LTD.','pharose remedies Ltd.',
                                                          'Naman Pharma Drugs','Signature Phytochemical Industries',
                                                          'Skohind Labs strategic alliance with Rousel Laboratories',
                                                          'uni medicolabs','Medreich limited','Centurion Laboratories PVT. LTD.'))
index_man <- fread("input/renumbering_manuf_index.txt", sep="\t", header=TRUE)
df <- merge(index_man, df, by.x="manufacturer", by.y="ATTRIBUTE_Info_Label_Manufacturer")



# Plots
plot_value <- df$ATTRIBUTE_API

score_plot <- ggplot(df, aes(as.numeric(PC1), as.numeric(PC3)))+
  geom_point(aes(colour=as.factor(plot_value)), pch=16, cex=1.0, alpha=0.7) +
  #geom_point(colour="black", pch=16, cex=1.0, alpha=0.7) +
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  #geom_text(aes(label=filename), size = 2)+
  theme_minimal()+
  theme(panel.background = element_rect(fill="grey85"),
        panel.grid.major=element_line(colour="grey90",size=0.25, linetype="dashed"),
        panel.grid.minor=element_line(colour="grey90",size=0.25, linetype="dashed"),
        panel.border = element_rect(colour ="black", size=0.25, linetype = "solid", fill = NA),
        axis.ticks=element_line(colour ="black",size=0.25, linetype="solid"),
        axis.text=element_text(colour="black",size=6),
        axis.text.x=element_text(colour="black",size=6),
        axis.line=element_line(colour="black",size=0.25, linetype="solid"),                  
        axis.title=element_text(colour="black",size=6),
        strip.text.x = element_text(colour="black",size=6),
        aspect.ratio=1,
        legend.position="none",
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size=6)) +
  labs(x=paste("PC1 (",round(mypca@R2[1] *100,2), "%)"), y=paste("PC3 (",round(mypca@R2[3] *100,2), "%)"))
print(score_plot)
ggsave(score_plot, filename="R_output/Drugs_FUSED_PCA_PC1vsPC3_Score_All_API.pdf", width = 2, height = 2, units = "in", dpi=300, useDingbats=FALSE)

# Plots
plot_value <- df$ATTRIBUTE_API_Drug2

score_plot <- ggplot(df, aes(as.numeric(PC1), as.numeric(PC3)))+
  geom_point(aes(colour=as.factor(plot_value)), pch=16, cex=1.0, alpha=0.7) +
  scale_color_brewer(palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  #geom_text(aes(label=filename), size = 2)+
  theme_minimal()+
  theme(panel.background = element_rect(fill="grey85"),
        panel.grid.major=element_line(colour="grey90",size=0.25, linetype="dashed"),
        panel.grid.minor=element_line(colour="grey90",size=0.25, linetype="dashed"),
        panel.border = element_rect(colour ="black", size=0.25, linetype = "solid", fill = NA),
        axis.ticks=element_line(colour ="black",size=0.25, linetype="solid"),
        axis.text=element_text(colour="black",size=6),
        axis.text.x=element_text(colour="black",size=6),
        axis.line=element_line(colour="black",size=0.25, linetype="solid"),                  
        axis.title=element_text(colour="black",size=6),
        strip.text.x = element_text(colour="black",size=6),
        aspect.ratio=1,
        legend.position="none",
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size=6)) +
  labs(x=paste("PC1 (",round(mypca@R2[1] *100,2), "%)"), y=paste("PC3 (",round(mypca@R2[3] *100,2), "%)"))
print(score_plot)
ggsave(score_plot, filename="R_output/Drugs_FUSED_PCA_PC1vsPC3_Score_All_APIDrug2.pdf", width = 2, height = 2, units = "in", dpi=300, useDingbats=FALSE)


# facet on Manufacturer and Package ID
score_plot <- ggplot(df, aes(as.numeric(PC1), as.numeric(PC3)))+
  geom_point(data = select(df, -c(number_man,Cbind_Sample_ID)), aes(colour=ATTRIBUTE_API), pch=16, cex=0.5, alpha=0.5) +
  geom_point(aes(fill=ATTRIBUTE_API), colour="black", pch=21, cex=1.0, alpha=1.0) +
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  theme_minimal()+
  facet_wrap(~number_man~Cbind_Sample_ID, ncol = 7)+
  theme(panel.background = element_rect(fill="grey85"),
        panel.grid.major=element_line(colour="grey90",size=0.25, linetype="dashed"),
        panel.grid.minor=element_line(colour="grey90",size=0.25, linetype="dashed"),
        panel.border = element_rect(colour ="black", size=0.25, linetype = "solid", fill = NA),
        axis.ticks=element_line(colour ="black",size=0.25, linetype="solid"),
        axis.text=element_text(colour="black",size=6),
        axis.text.x=element_text(colour="black",size=6),
        axis.line=element_line(colour="black",size=0.25, linetype="solid"),                  
        axis.title=element_text(colour="black",size=6),
        strip.text.x = element_text(colour="black",size=6),
        aspect.ratio=1,
        legend.position="none",
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size=6)) +
  labs(x=paste("PC1 (",round(mypca@R2[1] *100,2), "%)"), y=paste("PC3 (",round(mypca@R2[3] *100,2), "%)"))
print(score_plot)
ggsave(score_plot, filename="R_output/Drugs_FUSED_PCA_PC1vsPC3_Score_All_API_Facet.pdf", width = 7.1, height = 11, units = "in", dpi=300, useDingbats=FALSE)


# variance plots
variance <- round(mypca@R2*100,2)
cumul_variance <- round(mypca@R2cum*100,2)
PCs_num <- c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8")
plot_variance <- data.frame(PCs_num, variance, cumul_variance)

# facet on Manufacturer and Package ID
variance <- ggplot(plot_variance, aes(as.numeric(PCs_num), as.numeric(variance)))+
  geom_point(colour="black", pch=16, cex=1.25, alpha=1.0)+
  geom_path(aes(as.numeric(PCs_num), as.numeric(variance)),colour="black")+
  geom_point(aes(x=as.numeric(PCs_num), y=as.numeric(cumul_variance)), colour="red", pch=1, cex=1.25, alpha=1.0)+
  geom_path(aes(as.numeric(PCs_num), as.numeric(cumul_variance)),colour="red")+
  theme_minimal()+
  theme(panel.background = element_rect(fill="grey85"),
        panel.grid.major=element_line(colour="grey90",size=0.25, linetype="dashed"),
        panel.grid.minor=element_line(colour="grey90",size=0.25, linetype="dashed"),
        panel.border = element_rect(colour ="black", size=0.25, linetype = "solid", fill = NA),
        axis.ticks=element_line(colour ="black",size=0.25, linetype="solid"),
        axis.text=element_text(colour="black",size=6),
        axis.text.x=element_text(colour="black",size=6),
        axis.line=element_line(colour="black",size=0.25, linetype="solid"),                  
        axis.title=element_text(colour="black",size=6),
        strip.text.x = element_text(colour="black",size=6),
        aspect.ratio=1,
        legend.position="none",
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size=6)) +
  labs(x="PCs", y="% variance")
print(variance)
ggsave(variance, filename="R_output/Drugs_FUSED_PCA_Variance.pdf", width = 2.25, height = 2.25, units = "in", dpi=300, useDingbats=FALSE)




