library(naturalsort)
library("corrplot")
# library(ggrgl)
library(tidyverse)
library(ggsci)
library(pcaMethods)
library("FactoMineR")
library("factoextra")
library("PerformanceAnalytics")
library(ggpubr)
library(ggrepel)
library(ggforce)
library(sp)
library(fasano.franceschini.test)
library(patchwork)

ks_test_2d <- function(data, col_name, res_pca){
	groups <- names(table(data_meta[[col_name]]))
	pair_groups <- combn(groups, 2)
	ks_test_2d <- apply(pair_groups, 2, function(pair_t) {
		dim1_t <- res_pca$ind$coord[,1]
		dim2_t <- res_pca$ind$coord[,2]
		
		dim1_t1 <- dim1_t[data_meta[[col_name]] == pair_t[1]]
		dim1_t2 <- dim1_t[data_meta[[col_name]] == pair_t[2]]
		dim2_t1 <- dim2_t[data_meta[[col_name]] == pair_t[1]]
		dim2_t2 <- dim2_t[data_meta[[col_name]] == pair_t[2]]

		g1Data <- data.frame(x = dim1_t1, y = dim2_t1)
		g1Data <- g1Data[!is.na(g1Data[,1]),]
		g2Data <- data.frame(x = dim1_t2, y = dim2_t2)
		g2Data <- g2Data[!is.na(g2Data[,1]),]

		ks_test_2d <- fasano.franceschini.test(g1Data, g2Data)
		# str(ks_test_2d)
		ks_test_2d$p.value
	})
	df_ks_test <- as_tibble(as.data.frame(t(pair_groups)))
	df_ks_test$p_value <- ks_test_2d
	df_ks_test$`<0.05?` <- df_ks_test$p_value < 0.05
	df_ks_test$`<0.01?` <- df_ks_test$p_value < 0.01
	names(df_ks_test)[1:2] <- c("Group_1", "Group_2")
	df_ks_test$variable <- col_name
	return(df_ks_test)
}

wilcox_test_d1 <- function(data, col_name, res_pca){
	groups <- names(table(data_meta[[col_name]]))
	pair_groups <- combn(groups, 2)
	wilc_test_2d <- apply(pair_groups, 2, function(pair_t) {
		dim1_t <- res_pca$ind$coord[,1]
		
		dim1_t1 <- dim1_t[data_meta[[col_name]] == pair_t[1]]
		dim1_t2 <- dim1_t[data_meta[[col_name]] == pair_t[2]]
		
		rst_test <- wilcox.test(dim1_t1, dim1_t2)
		# str(ks_test_2d)
		rst_test$p.value
	})
	df_wilc_test <- as_tibble(as.data.frame(t(pair_groups)))
	df_wilc_test$p_value <- wilc_test_2d
	df_wilc_test$`<0.05?` <- df_wilc_test$p_value < 0.05
	df_wilc_test$`<0.01?` <- df_wilc_test$p_value < 0.01
	names(df_wilc_test)[1:2] <- c("Group_1", "Group_2")
	return(df_wilc_test)
}

data_raw <- read_csv("../data/data_raw_filled.csv")
data_raw$Age_group[data_raw$Age_group=="X"] <- "Non-infected children"
data_raw$Age_group[data_raw$Age_group=="Y"] <- "Non-infected adults"
data_m <- data_raw %>% select(IgG_S:`FcgRIIa Avidity Index_OC43`)
unique(apply(data_m, 2, function(x){which(is.na(x))})) # row 82 had missing value
data <- data_raw[-82,]

# inital PCA
data_meta <- data %>% select(-(IgG_S:`FcgRIIa Avidity Index_OC43`))
data_m <- data %>% select(`IgG_S`:`FcgRIIa Avidity Index_OC43`)

data_m <- data_m %>% mutate_all(as.numeric)
data_completed <- prep(data_m, scale= "uv")

## PCA
cba_cur="initial"
res.pca <- PCA(data_completed, scale.unit = TRUE, ncp = 6, graph = FALSE)
df <- as_tibble(res.pca$eig)
df <- bind_cols(PC=rownames(res.pca$eig), df)
(p0 <- fviz_screeplot(res.pca, ncp=14, addlabels = TRUE))
file_out <- paste0("../results/screet_plot_", cba_cur, ".pdf")
ggsave(file_out, width = 8, height = 6)

(p_tmp <- fviz_pca_var(res.pca, gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, select.var = list(contrib=15)))

(p <- fviz_contrib(res.pca, "var", axes = c(1,2)))
file_out <- paste0("../results/contributions_", cba_cur, ".pdf")
ggsave(file_out, width = 8, height = 5)

data_high_con <- p$data %>% filter(contrib>=1/nrow(p$data)*100)
data_scree <- p0$data
data_scree$label <- paste0(rownames(data_scree), " (", round(data_scree$eig,2), "%)")

df_plot <- as.data.frame(res.pca$var$coord)
df_plot$label <- rownames(df_plot)
df_plot$color <- ifelse(df_plot$label %in% as.character(data_high_con$name), "#ca0020", "#bababa")

df_plot$color[df_plot$label %in% c("FcgRIIIa_OC43", "IgG_OC43", "IgA_OC43", "FcgRIIa_OC43")] <- "#0571b0"
df_plot$color_seg <- ifelse(df_plot$color=="#bababa", "#bababa", "#000000")

ggplot(df_plot)+
	geom_circle(aes(x0=0, y0=0, r=1))+
	theme_minimal()+
	geom_hline(yintercept=0, linetype="dashed")+
	geom_vline(xintercept=0, linetype="dashed")+
	geom_segment(aes_string(x = 0, y = 0, xend = "Dim.1", yend = "Dim.2", color="color_seg"), arrow = grid::arrow(length = grid::unit(0.2,"cm")))+
	geom_text_repel(aes(x=Dim.1, y=Dim.2, label=label, color=color), size=2, bg.color = "white", segment.size=1, data=. %>% filter(color!="#bababa"))+
	geom_text_repel(aes(x=Dim.1, y=Dim.2, label=label), color="#bababa", size=2, bg.color = "white", segment.size=1, data=. %>% filter(color=="#bababa"))+
	scale_color_identity()+
	xlab(data_scree$label[1])+
	ylab(data_scree$label[2])+	
	# coord_fixed()+
	NULL
file_out <- paste0("../results/var_plot_", cba_cur, ".pdf")
ggsave(file_out, width = 8, height = 8)

## symptoms
data_meta$Symptoms <- ifelse(data_meta$Symptoms==1, "Symptomatic", "Asymptomatic")
data_meta$Symptoms[is.na(data_meta$Symptoms)] <- "Non-infected"
data_meta$Symptoms <- factor(data_meta$Symptoms, levels=c("Symptomatic", "Asymptomatic", "Non-infected"))
(p_symptoms <- fviz_pca_ind(res.pca, geom.ind = "point", col.ind = data_meta$Symptoms, palette = "uchicago", addEllipses = TRUE, legend.title = "Symptoms", title="Symptoms"))
file_out <- paste0("../results/idv_by_symptom_", cba_cur, ".pdf")
ggsave(file_out, width = 8, height = 6)
(df_ks_test_2d_symptoms <- ks_test_2d(data_meta, "Symptoms", res.pca))
write_csv(df_ks_test_2d_symptoms, paste0("../results/idv_by_symptom_", cba_cur, ".csv"))
wilcox_test_d1(data_meta, "Symptoms", res.pca)

## Age
(p_age <- fviz_pca_ind(res.pca, geom.ind = "point", col.ind = data_meta$Age_group, palette = "uchicago", addEllipses = TRUE, legend.title = "Age group", title="Age"))
file_out <- paste0("../results/idv_by_age_", cba_cur, ".pdf")
ggsave(file_out, width = 8, height = 6)
(df_ks_test_2d_agegroup <- ks_test_2d(data_meta, "Age_group", res.pca))
write_csv(df_ks_test_2d_agegroup, paste0("../results/idv_by_age_", cba_cur, ".csv"))
wilcox_test_d1(data_meta, "Age_group", res.pca)

## Sex
(p_sex <- fviz_pca_ind(res.pca, geom.ind = "point", col.ind = data_meta$Sex, palette = "uchicago", addEllipses = TRUE, legend.title = "Sex", title="Sex"))
file_out <- paste0("../results/idv_by_sex_", cba_cur, ".pdf")
ggsave(file_out, width = 8, height = 6)
(df_ks_test_2d_sex <- ks_test_2d(data_meta, "Sex", res.pca))
(df_ks_test_2d_sex <- ks_test_2d(data_meta, "Sex", res.pca))
write_csv(df_ks_test_2d_agegroup, paste0("../results/idv_by_sex_", cba_cur, ".csv"))
wilcox_test_d1(data_meta, "Sex", res.pca)

# ## Group
# (p_tmp <- fviz_pca_ind(res.pca, geom.ind = "point", col.ind = data_meta$Group, palette = "uchicago", addEllipses = TRUE, legend.title = "Sex"))
# (df_ks_test_2d_group <- ks_test_2d(data_meta, "Group", res.pca))
# wilcox_test_d1(data_meta, "Group", res.pca)

## Age symptoms
data_meta$Age_group_more <- paste0(data_meta$Age_group, "-", data_meta$Symptoms)
data_meta$Age_group_more <- gsub("-Non-infected", "", data_meta$Age_group_more)
(p_tmp <- fviz_pca_ind(res.pca, geom.ind = "point", col.ind = data_meta$Age_group_more, palette = "uchicago", addEllipses = TRUE, legend.title = "Age group"))
file_out <- paste0("../results/idv_by_age_more_", cba_cur, ".pdf")
ggsave(file_out, width = 8, height = 6)
(df_ks_test_2d_age_group_more <- ks_test_2d(data_meta, "Age_group_more", res.pca))
wilcox_test_d1(data_meta, "Age_group_more", res.pca)

## Time
(p_time <- fviz_pca_ind(res.pca, geom.ind = "point", col.ind = data_meta$Time, gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), legend.title = "Time post infection"))
data_meta$time_group <- ifelse(data_meta$Time<=14, "<=14 days", ">14 days")
data_meta$time_group[is.na(data_meta$time_group)] <- "Non-infected"
(p_time <- fviz_pca_ind(res.pca, geom.ind = "point", col.ind = data_meta$time_group, palette = "uchicago", addEllipses = TRUE, legend.title = "Acute infection", title="Time"))
file_out <- paste0("../results/idv_by_timegroup_", cba_cur, ".pdf")
ggsave(file_out, width = 8, height = 6)
(df_ks_test_2d_timegroup <- ks_test_2d(data_meta, "time_group", res.pca))
write_csv(df_ks_test_2d_timegroup, paste0("../results/idv_by_timegroup_", cba_cur, ".csv"))
wilcox_test_d1(data_meta, "time_group", res.pca)

(p_symptoms | p_age)/
(p_sex | p_time)
file_out <- paste0("../results/idv_by_all_", cba_cur, ".pdf")
ggsave(file_out, width = 16, height = 12)

df_test_all <- bind_rows(df_ks_test_2d_symptoms, df_ks_test_2d_agegroup, df_ks_test_2d_sex, df_ks_test_2d_timegroup)
df_test_all <- df_test_all %>% select(variable, everything())
write_csv(df_test_all, "../results/idc_by_all_test.csv")

qplot(data_meta$Time, res.pca$ind$coord[,1], color = Age_group, data=data_meta)
qplot(data_meta$Time, res.pca$ind$coord[,2], color = Age_group, data=data_meta)

## Age time
data_meta$Age_group_time <- paste0(data_meta$Age_group, " (", ifelse(data_meta$Time<=14, "<=14 days", ">14 days"), ")")
data_meta$Age_group_time <- gsub(" (NA)", "", data_meta$Age_group_time, fixed=T)
(p_tmp <- fviz_pca_ind(res.pca, geom.ind = "point", col.ind = data_meta$Age_group_time, palette = "uchicago", addEllipses = TRUE, legend.title = "Age group (Time)"))
file_out <- paste0("../results/idv_by_age_time_", cba_cur, ".pdf")
ggsave(file_out, width = 8, height = 6)


# # Second PCA
# cba_cur="second"
# data_meta <- data %>% select(-(IgG_S:`FcgRIIa Avidity Index_OC43`))
# data <- data %>% filter(!grepl("Non-infected", Age_group))
# data_m <- data %>% select(`IgG_S`:`FcgRIIa Avidity Index_OC43`)
# data_m <- data_m %>% select(!contains("ratio"))
# data_m <- data_m %>% select(!contains("Avidity"))
# data_m <- data_m %>% mutate_all(as.numeric)
# data_completed <- prep(data_m, scale= "uv")

# ## PCA
# res.pca <- PCA(data_completed, scale.unit = TRUE, ncp = 6, graph = FALSE)
# df <- as_tibble(res.pca$eig)
# df <- bind_cols(PC=rownames(res.pca$eig), df)
# fviz_screeplot(res.pca, ncp=14, addlabels = TRUE)
# file_out <- paste0("../results/screet_plot_", cba_cur, ".pdf")
# ggsave(file_out, width = 8, height = 6)

# (p_tmp <- fviz_pca_var(res.pca, gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, select.var = list(contrib=15)))

# (p <- fviz_contrib(res.pca, "var", axes = c(1,2)))
# file_out <- paste0("../results/contributions_", cba_cur, ".pdf")
# ggsave(file_out, width = 8, height = 5)

# ## symptoms
# data_meta$Symptoms <- ifelse(data_meta$Symptoms==1, "Symptomatic", "Asymptomatic")
# data_meta$Symptoms[grepl("Non-infected", data_meta$Age_group)] <- "Non-infected"
# (p_tmp <- fviz_pca_ind(res.pca, geom.ind = "point", col.ind = data_meta$Symptoms, palette = "uchicago", addEllipses = TRUE, legend.title = "Symptoms"))
# file_out <- paste0("../results/idv_by_symptom_", cba_cur, ".pdf")
# ggsave(file_out, width = 8, height = 6)

# ## Age
# (p_tmp <- fviz_pca_ind(res.pca, geom.ind = "point", col.ind = data_meta$Age_group, palette = "uchicago", addEllipses = TRUE, legend.title = "Age group"))
# file_out <- paste0("../results/idv_by_age_", cba_cur, ".pdf")
# ggsave(file_out, width = 8, height = 6)

# ## Sex
# (p_tmp <- fviz_pca_ind(res.pca, geom.ind = "point", col.ind = data_meta$Sex, palette = "uchicago", addEllipses = TRUE, legend.title = "Sex"))

# ## Group
# (p_tmp <- fviz_pca_ind(res.pca, geom.ind = "point", col.ind = data_meta$Group, palette = "uchicago", addEllipses = TRUE, legend.title = "Sex"))

# ## Age symptoms
# data_meta$Age_group_more <- paste0(data_meta$Age_group, "-", data_meta$Symptoms)
# data_meta$Age_group_more <- gsub("-Non-infected", "", data_meta$Age_group_more)
# (p_tmp <- fviz_pca_ind(res.pca, geom.ind = "point", col.ind = data_meta$Age_group_more, palette = "uchicago", addEllipses = TRUE, legend.title = "Age group"))
# file_out <- paste0("../results/idv_by_age_more_", cba_cur, ".pdf")
# ggsave(file_out, width = 8, height = 6)

# ## Time
# (p_tmp <- fviz_pca_ind(res.pca, geom.ind = "point", col.ind = data_meta$Time, gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), legend.title = "Time"))

# qplot(data_meta$Time, res.pca$ind$coord[,1], color = Age_group, data=data_meta)
# qplot(data_meta$Time, res.pca$ind$coord[,2], color = Age_group, data=data_meta)

# ## Age time
# data_meta$Age_group_time <- paste0(data_meta$Age_group, " (", ifelse(data_meta$Time<=14, "<=14", ">14"), ")")
# data_meta$Age_group_time <- gsub(" (NA)", "", data_meta$Age_group_time, fixed=T)
# (p_tmp <- fviz_pca_ind(res.pca, geom.ind = "point", col.ind = data_meta$Age_group_time, palette = "uchicago", addEllipses = TRUE, legend.title = "Age group"))
# file_out <- paste0("../results/idv_by_age_time_", cba_cur, ".pdf")
# ggsave(file_out, width = 8, height = 6)






# test cytokine diff between severe and normal
if(!any(is.infinite(cytokine_res))){
	ks_tst <- kruskal.test(cytokine_res, data_meta$Symptoms)
	tmp <- tapply(cytokine_res, data_meta$Symptoms, median)
	pair_t <- combn(names(tmp), 2)

	wc_rst_p <- sapply(seq_len(ncol(pair_t)), function(j){
		wc_rst <- wilcox.test(cytokine_res[data_meta$Symptoms == pair_t[1,j]], cytokine_res[data_meta$Symptoms == pair_t[2,j]])
		wc_rst$p.value	
	})
	df_ks_test <<- bind_rows(df_ks_test, tibble(cytokine_res = cba_cur, p_value_kruskal_test = ks_tst$p.value, median_asymp = tmp[1], median_severe = tmp[2], median_symp = tmp[3], p_value_wc_test_asym_seve = wc_rst_p[1], p_value_wc_test_asym_symp = wc_rst_p[2], p_value_wc_test_seve_symp = wc_rst_p[3]))
}
