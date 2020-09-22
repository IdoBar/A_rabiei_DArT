# load custom functions from github
devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "util.R")

# install dependencies
# BiocManager::install(c("SNPRelate", "qvalue"))
# install.packages(c("backports", "colorspace", "adegenet"))
# load packages we're going to use
# BiocManager::install("qvalue", update = FALSE)
pacman::p_load(  tidyverse, forcats, RColorBrewer, paletteer, ggrepel, glue, dartR, poppr, SNPassoc, here,
                 dendextend, ComplexHeatmap, pheatmap) # quadprog, raster, mvtnorm,pheatmap

#### Read DArT Data ####
pheno_file <- recent_file("data", "A_rabiei_pathotypes.+.xlsx")
isolate_table <- readxl::read_excel(path = pheno_file , sheet = "pathotyping") %>% #select(-one_of(c("Location"))) %>% 
  mutate(Taxa=Isolate) %>%
  mutate(Host=sub("^Genesis09[0-9]+$", "Genesis090", 
                  sub(" ", "", sub("PBA ", "", Host))))
taxa_regions <- readxl::read_excel("output/dart_samples_sum_IB_S.xlsx", sheet = "taxa_regions") %>% 
  select(Taxa, Region)

taxa_table <- readxl::read_excel("data/Taxa_DArTseq_Plates1-2_5-05-19.xlsx", na = c("N/A", "")) %>%
  inner_join(taxa_regions) %>% 
  mutate(Host=sub("^Genesis09[0-9]+$", "Genesis090", 
                  sub(" ", "", sub("PBA ", "", Host)))) %>% 
  mutate_at(.vars=c("Year"), as.numeric) %>% select(Taxa, Location,  State, Region, Year, Host, Haplotype) %>% 
  full_join(isolate_table, by = c("Taxa", "State", "Year", "Host")) %>% mutate(State=sub("[^A-Z]", "", toupper(State))) %>% 
  mutate(Location=if_else(is.na(Location.x), Location.y, Location.x)) 
  # Host=sub("CICA.+", "CICA1152", Host), 
write_xlsx(taxa_table, "data/5_year_complete_SSR_db.xlsx", "2013-2018_pheno", 
           overwritesheet=TRUE, asTable=FALSE)
  
taxa_table %>% filter(is.na(Location))
ssr_haplos <- readxl::read_excel("data/5_year_complete_SSR_db.xlsx")


# gl1_2 <- gl.read.dart("data/Report_DAsc19-4108_ArabieiPlates1-2/Report_DAsc19-4108_1_moreOrders_SNP_mapping_2.csv") 
gl3 <- gl.read.dart("data/Report-DAsc19-4353_ArabieiPlate3/Report_DAsc19-4353_1_moreOrders_SNP_mapping_2.csv")

# fix replaced names
replace_names <- c("122/17-1"= "122-1/17","122/17-2"="122-2/17","122/17-3"="122-3/17","14HOR15"="14HOR015",
                             "17HOR009"="17HOR002",
                   # "F14038"="14038", "F14039"="14039", "F14040"="14040", "F14090"="FT14090", "F17201-1"="17201-1",
                             "F15009"="FT15009","F15019"="FT15019", "F15021"="FT15021", 
                   "F16181-1"="FT16187-1", 
                             "F16188-1"="FT16188-1", 
                   # "MER-17-373"="MER17373", "MER-17-374"= "MER17374", "MER-17-376"= "MER17376", "MER-17-378"= "MER17378", "MER-17-379"= "MER17379", "MER-17-380"="MER17380", "MER-17-381"="MER17381",                 "MER-17-382"="MER17382", 
                   "AS17083"="AS18083", "AS17087"="AS18087")


taxa_table %>% count(Host) %>% filter(grepl("CICA", Host))
taxa_table %>% count(Taxa) %>% filter(n>1)
taxa_table %>% filter(grepl("MER", Taxa))
taxa_table %>% count(Region)

metadata <- tibble(indNames=indNames(gl3)) %>% 
  mutate(Taxa=if_else(is.na(replace_names[indNames]), indNames, replace_names[indNames]))  %>% 
  # left_join(isolate_table, by = c("Taxa"="Isolate")) %>% 
  left_join(taxa_table) 


# %>% 
#  mutate(Region=case_when(State=="QLD"~"Reg1", State=="WA"~"Reg6",State=="SA"~"Reg5", State=="VIC"~"Reg4", TRUE~Region))

metadata$Region[metadata$State=="WA"] <- "Reg6"
metadata %>% filter(is.na(Region))
  

metadata <- metadata %>% 
  select(id=indNames, Location, State, Region, Year, Host, Pathotype, Path_rating, Haplotype)  %>% 
  write_csv("data/DArT_metadata.csv")

# metadata <- read_csv(here("data/DArT_metadata.csv"))
# double check matching names and length
length(indNames(gl3)) == nrow(metadata)
all(indNames(gl3) == metadata$indNames)

# check how many isolate names don't appear in the pathotype table
metadata %>% filter(!Taxa %in% isolate_table$Taxa) %>% .[c("indNames", "Taxa")]

# check which isolates don't have haplotype information (from each year)
metadata %>% filter(is.na(Haplotype)) %>% count(Year)

# read in data with metadata in place
Arab_gl <- gl.read.dart("data/Report-DAsc19-4353_ArabieiPlate3/Report_DAsc19-4353_1_moreOrders_SNP_mapping_2.csv", 
                        ind.metafile = "data/DArT_metadata.csv")

sum_table <- read_csv("data/DArT_metadata.csv") %>% filter(State!="SPAIN") %>% count(Year, State) %>% 
  pivot_wider(names_from = State, values_from = n) %>% mutate_all(~replace_na(., 0)) %>% 
  mutate(Total=rowSums(.[-1]), Year=as.character(Year)) %>% 
  rbind(., c(NA, as.double(summarise_if(., is_double, ~sum(.))))) %>% mutate(Year = replace_na(Year, "Total"))

# set consistent strata
strata_factors <- Arab_gl@other$ind.metrics %>%
  # select(indNames, Location, State, Region, Year, Host, Pathotype, Path_rating, Haplotype, Taxa) %>%
  mutate(rowname=id, State=factor(State, levels = c("QLD", "NSW", "VIC", "SA", "WA", "SPAIN")),
         Year=factor(Year, levels=as.character(min(Year, na.rm = TRUE):max(Year, na.rm = TRUE))),
         Host=fct_relevel(fct_infreq(fct_lump_min(Host, min=6)), "Other", after = Inf),
         Haplotype=fct_relevel(fct_infreq(fct_lump_min(Haplotype, min=2)), "Other", after = Inf),
         Patho.Group=factor(Pathotype, levels=paste0("Group", 0:5))) %>%
  mutate_at(vars(Year, Host, Haplotype, Patho.Group, Region), ~fct_explicit_na(., na_level = "Unknown")) %>% 
  column_to_rownames()


strata_factors %>% count(Host)
#### Exploratory Data Analysis ####
ploidy(Arab_gl) <- 1
gl.report.monomorphs(Arab_gl)  # all loci are polymorphic
# check missing rates
gl.report.callrate(Arab_gl, method = "loc", v=3)
gl.report.callrate(Arab_gl, method = "ind", v = 3) # the 2 individuals from Spain are the most remote
# check reproducibility (RepAVG)
gl.report.repavg(Arab_gl)
# check multiple loci per tag and reproducibility
gl.report.secondaries(Arab_gl)

Arab_gl@strata <- strata_factors
setPop(Arab_gl) <- ~State
# apply filters for export
# find loci with high number of SNPs
snps_in_tag_thres <- 5
drop_loci <- Arab_gl@other$loc.metrics %>% count(CloneID) %>% filter(n>=snps_in_tag_thres)
# keep_loci <- Arab_gl@other$loc.metrics %>% count(CloneID) %>% filter(n<snps_in_tag_thres)

drop_snps <- locNames(Arab_gl)[sub("^(\\d+)-.+", "\\1", locNames(Arab_gl)) %in% drop_loci$CloneID]
# remove Spain samples and loci with more than 5 SNPs
# Arab_filt <- Arab_gl[!grepl("spain", indNames(Arab_gl), ignore.case = TRUE),
#                       sub("^(\\d+)-.+", "\\1", locNames(Arab_gl)) %in% keep_loci$CloneID, treatOther=TRUE] %>% 
#             gl.filter.monomorphs(v = 3) 
Arab_filt <- Arab_gl %>% gl.drop.loc(drop_snps, v=3) %>% 
  gl.drop.ind(indNames(Arab_gl)[grepl("spain", indNames(Arab_gl), ignore.case = TRUE)], v=3) 
  
dim(Arab_filt@other$loc.metrics)
gl.report.monomorphs(Arab_filt)
# Arab_filt@other$loc.metrics <- Arab_filt@other$loc.metrics %>% filter(CloneID %in% keep_loci$CloneID)
missing_thres <- 0.3

glsub <- Arab_filt %>% gl.filter.callrate(., method = "loc", threshold = (1-missing_thres), v = 3) %>% 
  gl.filter.callrate(., method = "ind", threshold = (1-missing_thres), recalc = TRUE, v = 3) %>% 
  gl.filter.repavg(., threshold = 0.95, v = 3)

ploidy(glsub) <- 1
glsub@strata <- strata_factors %>% filter(id %in% indNames(glsub)) %>% mutate_if(is.factor, ~fct_drop(.))

# calculate AMOVA (with clone correction)
agc <- as.genclone(gl2gi(glsub))
agc@strata <- strata_factors %>% filter(id %in% indNames(agc)) %>% mutate_if(is.factor, ~fct_drop(.))
amova.result <- poppr.amova(agc, ~State/Year, clonecorrect = TRUE)
amova.cc.test <- randtest(amova.result)
plot(amova.cc.test)
amova.cc.test


# define set colour palletes for each factor
# Year
# Host
# State
# Pathogenicity

# Calculate Tree, heatmap and clusters ####
tree <- aboot(glsub, tree = "upgma", distance = "bitwise.dist", sample = 100, showtree = T, 
              cutoff = 50, quiet = T,
              missing="ignore")
dist_matrix <- bitwise.dist(glsub, euclidean = TRUE)
# define clusters
my_hclust <- hclust(dist_matrix, method = "average")
cluster_num <- 6
mat_clusters <- cutree(tree = as.dendrogram(my_hclust), k = cluster_num) %>% tibble(id=names(.), Cluster=.)
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

# data_norm <- t(apply(dist_matrix, 1, cal_z_score))
# specify annotation for columns and rows
row_annotation <- glsub@strata  %>% select(id, State, Year) %>%  column_to_rownames("id") # mutate_all(~forcats::fct_drop(.)) %>%
col_annotation <- glsub@strata  %>% select(id, Patho.Group) %>% column_to_rownames("id") # inner_join(mat_clusters) %>% 

pc_colors <-  c("grey", brewer.pal(9, "RdYlGn")[round(seq(1,8, length.out = 9-1 ) ,0)]) %>% rev(.) %>%
  setNames(popNames(Arab_gclone))
# specify colours
paletteer_d("RColorBrewer::RdYlGn", direction = -1)
my_colours = list(
  State = levels(row_annotation$State) %>% setNames(as.character(paletteer_d("RColorBrewer::Set1", length(.))), .),
  Year = levels(row_annotation$Year) %>% setNames(as.character(paletteer_d("ggsci::default_uchicago", length(.))), .),
  Patho.Group = levels(col_annotation$Patho.Group) %>% 
    setNames(c(as.character(paletteer_d("RColorBrewer::RdYlGn", direction = -1))[round(seq(4,11, length.out = length(.)-1 ) ,0)], "grey80"), .)
)

# define colour function fro ComplexHeatmap
# col_fun <- circlize::colorRamp2(0:3, paletteer_c("viridis::viridis", 4, direction = -1))
# Heatmap(as.matrix(dist_matrix), name = "Genetic Distance", show_row_names = FALSE, show_column_names = FALSE, 
#         col = col_fun)


dist_heatmap <- pheatmap(as.matrix(dist_matrix), show_rownames = FALSE, show_colnames = FALSE, 
         color = paletteer_c("viridis::viridis", 50, direction = -1), annotation_colors = my_colours,
        # clustering_distance_rows = dist_matrix, clustering_distance_cols = dist_matrix,
         annotation_row = row_annotation, annotation_col = col_annotation, 
         cutree_rows = cluster_num, cutree_cols = cluster_num , 
        filename = filedate(glue("A_rabiei_samples_heatmap_{cluster_num}clusters"), ext=".pdf", outdir = here::here("output/plots")), width = 8, height = 6.5)

# calculate "enrichment" in clusters
high_path_thresh <- 3
heatmap_clusters <- cutree(tree = as.dendrogram(dist_heatmap$tree_row), k = cluster_num) %>% tibble(id=names(.), Cluster=.) 
pathog_cluster <- glsub@strata  %>% select(id, Path_rating) %>% 
  mutate(Virulence=case_when(Path_rating<high_path_thresh~"Low", Path_rating>=high_path_thresh~"High", TRUE~"Unknown")) %>% filter(Virulence!="Unknown") %>% 
  inner_join(heatmap_clusters) %>% mutate(Cluster=factor(LETTERS[Cluster]))

glsub@other$ind.metrics

# identify clusters on heatmap
clust_annotation <- glsub@strata  %>% select(id, Patho.Group)  %>% inner_join(heatmap_clusters) %>% 
  mutate(Cluster=factor(LETTERS[Cluster])) %>% column_to_rownames("id") 
# define colours
heatmap_pal <- my_colours
heatmap_pal$Cluster <- levels(clust_annotation$Cluster) %>% setNames(as.character(paletteer_d("rcartocolor::Bold", length(.))), .) # "rcartocolor::Prism"

# palettes_d_names %>% filter(type=="qualitative", length>=6, length<10)

pheatmap(as.matrix(dist_matrix), show_rownames = FALSE, show_colnames = FALSE, 
         color = paletteer_c("viridis::viridis", 50, direction = -1), annotation_colors = heatmap_pal,
         # clustering_distance_rows = dist_matrix, clustering_distance_cols = dist_matrix,
         annotation_row = row_annotation, annotation_col = clust_annotation , 
         cutree_rows = cluster_num, cutree_cols = cluster_num, 
         filename = filedate(glue("A_rabiei_samples_heatmap_{cluster_num}cluster_cols"), ext=".pdf",
         outdir = here("output/plots")),
width = 8, height = 6.5)
# make dataset to apply GLM
pacman::p_load(lme4, lsmeans, ggpubr, ggthemr, gtools)

# contingency table
cont_table <- table(pathog_cluster$Virulence, pathog_cluster$Cluster)
prop_table <- prop.table(table(pathog_cluster$Cluster, pathog_cluster$Virulence), margin = 1) %>% 
  as.data.frame() %>% rename(Cluster=Var1, pathog_class=Var2) %>% as_tibble() %>% 
  mutate(pathog_class=factor(pathog_class, levels = c("Low", "High")))


# save table
prop_table %>% pivot_wider(names_from = pathog_class, values_from = Freq) %>%  
  mutate_if(is.double, ~scales::percent(., accuracy = .1)) %>% 
  write_xlsx(., "output/A_rabiei_DArT_poppr.xlsx", sheet = "pathog_clusters2", asTable=FALSE, overwritesheet = TRUE)
# chisq.test(table(patho_cluster$Virulence, patho_cluster$Cluster))


# heatmap_clusters %>% mutate(value=1) %>% pivot_wider(names_from = Cluster, values_from = value, values_fill = list(value=0)) %>% 
#   pivot_longer()
# Cluster statistical analysis #####
cluster_data <- pathog_cluster %>% mutate(Cluster=LETTERS[Cluster], Virulence=ifelse(Path_rating<high_path_thresh,0, 1)) %>% write_csv(here("data/cluster_data.csv"))
data.glmer <- glmer(Virulence ~ Cluster + (1|id), family = binomial, data = cluster_data)
lsm <- lsmeans(data.glmer, "Cluster", type = "response")
signif_lsm <- pairs(lsm) %>% as.data.frame() %>% mutate(stars=stars.pval(p.value)) %>% 
  write_xlsx(., here("output/A_rabiei_DArT_poppr.xlsx"), sheet = "pathogenicity_cluster_OR_stats", asTable=FALSE, overwritesheet = TRUE) #%>% filter(p.value<0.05) 
lsm_mat <- matrix(NA, nrow = 6, ncol = 6,  dimnames = list(LETTERS[1:6], LETTERS[1:6])) 
lsm_mat[lower.tri(lsm_mat)] <- signif_lsm$odds.ratio
pairs(regrid(lsm)) %>% as.data.frame() %>% mutate(stars=stars.pval(p.value)) %>% 
  write_xlsx(., "output/A_rabiei_DArT_poppr.xlsx", sheet = "pathogenicity_cluster_diff_stats", asTable=FALSE, overwritesheet = TRUE)
# plot 
ggthemr("pale", text_size = 18)
# Stacked bar plots, add, based on paired ddelta
ggplot(prop_table, aes(x=Cluster, y=Freq, fill=pathog_class)) +
  geom_bar(stat="identity", width=0.6) +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  labs(fill="Virulence\nClassification") +
  annotate("text", x=levels(prop_table$Cluster), y = 1.05, label=c("ab", "ab", "c", "bc", "a", "ab"))
  # annotate("text", x=levels(prop_table$Cluster), y = 1.05, label=c("ab", "ab", "c\n***", "bc\n*", "a", "ab"))
ggsave(here::here("output/plots/A_rabiei_cluster_path_classification_bar.pdf"), width = 7, height=5)

# Stacked bar plots, based on odds-ratio
ggplot(prop_table, aes(x=Cluster, y=Freq, fill=pathog_class)) +
  geom_bar(stat="identity", width=0.6) +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  labs(fill="Virulence\nClassification") +
  geom_text(data = prop_table %>% filter(pathog_class == "High"), aes(label=scales::percent(Freq, accuracy = 3)), nudge_y = -0.05, colour="white", size=5) +
  annotate("text", x=levels(prop_table$Cluster), y = 1.05, label=c("ab", "ab", "c", "bc", "a", "ab"), size=5)
  # annotate("text", x="C", y = 1.05, label='***')
ggsave("output/plots/A_rabiei_cluster_path_classification_OR_signif.pdf", width = 7, height=5)

# visualise contingency table
# ggballoonplot(as.data.frame(cont_table) , fill = "value", color = "lightgray",
#               size = 10, show.label = TRUE) +
#   gradient_fill(c("blue", "white", "red"))
# add strata
# strata_factors <- metadata %>% filter(indNames %in% indNames(glsub)) %>%
#   select(indNames, Location, State, Region, Year, Host, Pathotype, Path_rating, Haplotype, Taxa) %>%
#   mutate(State=factor(State, levels = c("QLD", "NSW", "VIC", "SA", "WA")),
#          Year=factor(Year, levels=as.character(min(Year, na.rm = TRUE):max(Year, na.rm = TRUE))),
#          Host=fct_relevel(fct_infreq(fct_lump_min(Host, min=6)), "Other", after = Inf),
#          Haplotype=fct_relevel(fct_infreq(fct_lump_min(Haplotype, min=2)), "Other", after = Inf),
#          Pathotype=factor(Pathotype, levels=paste0("Group", 0:5))) %>%
#   mutate_at(vars(Year, Host, Haplotype, Pathotype, Region), ~fct_explicit_na(., na_level = "Unknown"))
# glsub@strata <- metadata %>% filter(indNames %in% indNames(glsub)) %>% select(Location, State, Region, Year, Host, Pathotype, Path_rating, Haplotype)
# setPop(glsub) <- ~Region
# export to faststructure
gl2faststructure(glsub, outfile = "output/fs/A_rabiei_fastructure.txt")
# export to findradstructure
head(glsub@other$loc.metrics)
# gl2fineradstructure <- function(gl){
#   # summarize metadata per tag (locus)
#   # find number of SNPs per tag (locus), collapse them with ;
#   # combine all allele combinations (haplotypes) - count unique and assign to each individual
#   
# }
tidy_gl <- radiator::genomic_converter(glsub, output = c("tidy"))
radiator::write_fineradstructure(tidy_gl$tidy.data, strata = glsub$strata %>% select(INDIVIDUALS = id, STRATA=State), 
                                 filename = "output/fs/A_rabiei_DArT")
# Make tree with all samples ####
# # apply filters for analysis
# missing_thres <- 0.4
# 
# glsub <- Arab_gl %>% gl.filter.callrate(., method = "loc", threshold = (1-missing_thres), v = 3) %>% 
#   gl.filter.callrate(., method = "ind", threshold = (1-missing_thres), recalc = TRUE, v = 3) %>% 
#   gl.filter.repavg(., threshold = 0.95, v = 3) %>% 
#   gl.filter.secondaries(., method = "best", v = 3)
#  
# ploidy(glsub) <- 1
# glsub@strata <- strata_factors %>% filter(id %in% indNames(glsub))
# # convert to genclone object
# # Arab_gclone <- as.genclone(gl2gi(glsub))
# # ploidy(Arab_gclone) <- 1
# # Arab_gclone@strata <- strata_factors %>% filter(id %in% indNames(Arab_gclone))
# # inherits(Arab_gclone, c("genlight", "genclone", "genind", "snpclone"))
# # create a tree
# tree <- aboot(glsub, tree = "upgma", distance = "bitwise.dist", sample = 100, showtree = T, cutoff = 50, quiet = T,
#               missing="ignore")
# dist_matrix <- bitwise.dist(glsub, euclidean = TRUE, mat=TRUE)
# 
# pheatmap(dist_matrix)




# apply filters for analysis of Australian Samples
# missing_thres <- 0.4
# gl_Aus <- Arab_gl[!grepl("spain", indNames(Arab_gl), ignore.case = TRUE),] %>% 
#   gl.filter.monomorphs(v = 3) %>% 
#   gl.filter.callrate(., method = "loc", threshold = (1-missing_thres), v = 3) %>% 
#   gl.filter.callrate(., method = "ind", threshold = (1-missing_thres), recalc = TRUE, v = 3) %>% 
#   gl.filter.repavg(., threshold = 0.95, v = 3) %>% 
#   gl.filter.secondaries(., method = "best", v = 3)
# 
# ploidy(gl_Aus) <- 1
# gl_Aus@strata <- strata_factors %>% filter(id %in% indNames(gl_Aus)) %>% mutate_if(is.factor, ~fct_drop(.))
# # Convert to allele frequencies (manage missing data)
# X <- adegenet::tab(gl_Aus, freq = TRUE, NA.method = "mean")
# pca1 <- ade4::dudi.pca(X, scale = FALSE, scannf = FALSE, nf=10)
# 
# pca_data <- pca1$li %>% rownames_to_column("id")  %>%
#   left_join(strata_factors) 
# # %>% mutate(Host=fct_relevel(fct_infreq(fct_lump_min(Host, min=6)), "Other", 
# #                                                   after = Inf),
# #                                  State=factor(State, levels = c("QLD", "NSW", "VIC", "SA", "WA", "SPAIN")))
# # Calculate variance for each component (columns starting with 'C')
# pca_var <- pca1$li %>% summarise_at(vars(starts_with("Axis")), var)
# # Calculate the percentage of each component of the total variance
# percentVar <- pca_var/sum(pca_var)
# # Define colour palette, but get rid of the awful yellow - number 6
# pal <- brewer.pal(9, "Set1")
# pal2 <- brewer.pal(8, "Dark2")
# text_cols <- adjustcolor( c("grey15", "dodgerblue3"), alpha.f = 0.8)
# # seq_face <- setNames(c("bold", "plain"),unique(pca_data$Sequenced))
# shapes <- c(21:25, 3)
# # Isolate to annotate
# outliers <- pca_data %>% filter(State=="Spain") 
# 
# #### PCA plot 1-2 geographic ####
# # shape - Year, Seq.Provider, Host (State)
# # Colors - Haplotype (qualitative), State (qualitative), Pathotype (diverging - RdYlGn)
# # Size - Pathogenicity
# # Create the plot for C1 and C2
# ggplot(pca_data, aes(x=Axis1, y=Axis2, fill=State)) + 
#   geom_point(alpha = 0.85, size=5, shape=21) +
#   scale_fill_paletteer_d(ggsci, category10_d3 ) +
#   plot_theme(baseSize = 20) + #  size="Year",
#   labs( x=glue::glue("C1: {round(percentVar[1]*100,2)}% variance"),
#         y=glue::glue("C2: {round(percentVar[2]*100,2)}% variance"))
# # Save plot to a pdf file
# ggsave(filedate(glue::glue("EDA_all_samples_PCA_state_PC1-2"),
#                 ext = ".pdf", outdir = "output/plots",
#                 dateformat = FALSE), width = 10, height=8)


#### Australian Pops Analysis ####
# filter on missing data (will remove Spanish isolates)
# apply filters
# missing_thres <- 0.1
# glsub <- Arab_gl %>% gl.filter.callrate(., method = "loc", threshold = (1-missing_thres), v = 3) %>% 
#   gl.filter.callrate(., method = "ind", threshold = (1-missing_thres), recalc = TRUE, v = 3) %>% 
#   gl.filter.repavg(., threshold = 0.95, v = 3) %>% 
#   gl.filter.secondaries(., method = "best", v = 3)

# Convert to allele frequencies (manage missing data)
X <- adegenet::tab(glsub, freq = TRUE, NA.method = "mean")
pca1 <- ade4::dudi.pca(X, scale = FALSE, scannf = FALSE, nf=10)

pca_data <- pca1$li %>% rownames_to_column("id")  %>%
  left_join(strata_factors)
# %>% mutate(Host=fct_relevel(fct_infreq(fct_lump_min(Host, min=6)), "Other", after = Inf),
#                                  State=factor(State, levels = c("QLD", "NSW", "VIC", "SA", "WA")))
# Calculate variance for each component (columns starting with 'C')
pca_var <- pca1$li %>% summarise_at(vars(starts_with("Axis")), var)
# Calculate the percentage of each component of the total variance
percentVar <- pca_var/sum(pca_var)
# Define colour palette, but get rid of the awful yellow - number 6
pal <- brewer.pal(9, "Set1")
pal2 <- brewer.pal(8, "Dark2")
text_cols <- adjustcolor( c("grey15", "dodgerblue3"), alpha.f = 0.8)
# seq_face <- setNames(c("bold", "plain"),unique(pca_data$Sequenced))
shapes <- c(21:25)
# Isolate to annotate
# replicated_samples <- pca_data %>% count(Isolate) %>% filter(n>1)
# outliers <- pca_data %>% filter(Axis2< -2.5) %>% .[,"Isolate"]
# outliers <- pca_data %>% filter(State=="Spain") 
# outliers <- c("16CUR017", "TR9543", "TR9529", "15CUR003") # PacBio-sequenced
# outliers <- pca_data$Isolate
#### PCA plots biological ####
# shape - Year, Seq.Provider, Host (State)
# Colors - Haplotype (qualitative), State (qualitative), Pathotype (diverging - RdYlGn)
# Size - Pathogenicity
# Create the plot for C1 and C2
plot_comps <- 1:2
plot_axes <- paste0("Axis", plot_comps)
ggplot(pca_data, aes(x=!!sym(plot_axes[1]), y=!!sym(plot_axes[2]), shape=State, fill=Patho.Group)) + 
  geom_point(alpha = 0.85, size=5, color="grey15") +
  scale_fill_manual(values = my_colours$Patho.Group) + 
  # scale_fill_brewer(palette = "RdYlGn", direction = -1) +
  # scale_fill_paletteer_d("RColorBrewer", "RdYlGn", direction = -1 ) + # ggsci, category10_d3
  scale_shape_manual(values = shapes) +
  # scale_color_manual(values=seq_cols) +
# # scale_size_continuous(range = c(3,6)) +
guides(shape = guide_legend(override.aes = list(size = 5 , fill=pal[2], color="grey15"), order = 1),  #, #
       fill = guide_legend(override.aes = list(shape = 21))) +   #  guide_legend(override.aes = list(size = 5, pch=shapes[1]), order = 2)
plot_theme(baseSize = 20) + #  size="Year",
  labs( x=glue::glue("C{plot_comps[1]}: {round(percentVar[plot_comps[1]]*100,2)}% variance"),
        y=glue::glue("C{plot_comps[2]}: {round(percentVar[plot_comps[2]]*100,2)}% variance"))
# Save plot to a pdf file
ggsave(filedate(glue::glue("Aus_samples_PCA_state_patho_PC{paste(plot_comps, collapse='-')}"),
                ext = ".pdf", outdir = here("output/plots")), width = 10, height=8)


# Create the plot for C3 and C4
plot_comps <- 3:4
plot_axes <- paste0("Axis", plot_comps)
ggplot(pca_data, aes(x=!!sym(plot_axes[1]), y=!!sym(plot_axes[2]), shape=State, fill=Patho.Group)) + 
  geom_point(alpha = 0.85, size=5, color="grey15") +
  scale_fill_manual(values = my_colours$Patho.Group) +  
  scale_shape_manual(values = shapes) +
  # scale_color_manual(values=seq_cols) +
  # # scale_size_continuous(range = c(3,6)) +
  guides(shape = guide_legend(override.aes = list(size = 5 , fill=pal[2], color="grey15"), order = 1),  #, #
         fill = guide_legend(override.aes = list(shape = 21))) +   #  guide_legend(override.aes = list(size = 5, pch=shapes[1]), order = 2)

  plot_theme(baseSize = 20) + #  size="Year",
  labs( x=glue::glue("C{plot_comps[1]}: {round(percentVar[plot_comps[1]]*100,2)}% variance"),
        y=glue::glue("C{plot_comps[2]}: {round(percentVar[plot_comps[2]]*100,2)}% variance"))
# Save plot to a pdf file
ggsave(filedate(glue::glue("Aus_samples_PCA_state_patho_PC{paste(plot_comps, collapse='-')}"),
                ext = ".pdf", outdir = here("output/plots")), width = 10, height=8)


#### Highlight high pathogenicity isolates  ####
# Create the plot for C1 and C2
plot_comps <- 3:4
plot_axes <- paste0("Axis", plot_comps)
path_cols <- c(rep(NA, 4),"#FDAE61" ,"#D73027", NA)

ggplot(pca_data, aes(x=!!sym(plot_axes[1]), y=!!sym(plot_axes[2]), shape=State, fill=Patho.Group)) + 
  geom_point(alpha = 0.85, size=5) +
  # scale_fill_brewer(palette = "RdYlGn", direction = -1) + 
  scale_fill_manual(values=path_cols) +
  scale_shape_manual(values = shapes) +
  # scale_color_manual(values=seq_cols) +
  # # scale_size_continuous(range = c(3,6)) +
  guides(shape = guide_legend(override.aes = list(size = 5 , fill=pal[2], color="grey15"), order = 1),  #, #
         fill = guide_legend(override.aes = list(shape = 21))) +   #  guide_legend(override.aes = list(size = 5, pch=shapes[1]), order = 2)
  # scale_size_manual(values = path_sizes[pca_data$Pathogenicity]) +
  
  # geom_text_repel(aes(label=Haplotype), size=4, point.padding = 0.75) +
  plot_theme(baseSize = 20) + #  size="Year",
  labs( x=glue::glue("C{plot_comps[1]}: {round(percentVar[plot_comps[1]]*100,2)}% variance"),
        y=glue::glue("C{plot_comps[2]}: {round(percentVar[plot_comps[2]]*100,2)}% variance"))
# Save plot to a pdf file
ggsave(filedate(glue::glue("Aus_samples_PCA_state_highpatho_PC{paste(plot_comps, collapse='-')}"),
                ext = ".pdf", outdir = here("output/plots")), width = 10, height=8)

#### PCA plots spatio-temporal ####
# shape - Year, Seq.Provider, Host (State)
# Colors - Haplotype (qualitative), State (qualitative), Pathotype (diverging - RdYlGn)
palettes_d_names %>% filter(type=="qualitative", length>11)
palettes_d_names %>% filter(type=="qualitative", length>11)
# Size - Pathogenicity
# Create the plot for C1 and C2
shapes <- 0:6# as.character(3:8)#
plot_comps <- 1:2
plot_axes <- paste0("Axis", plot_comps)
ggplot(pca_data, aes(x=!!sym(plot_axes[1]), y=!!sym(plot_axes[2]), shape=fct_inseq(factor(Year)), colour=State)) + 
  geom_point(alpha = 0.85, size=3, stroke = 1.5) +
  scale_colour_brewer(palette = "Set1") +
  # scale_colour_paletteer_d("rcartocolor", "Bold") + # ggsci, category10_d3
  scale_shape_manual(values = shapes) +
  # scale_color_manual(values=seq_cols) +
  # # scale_size_continuous(range = c(3,6)) +
  guides(shape = guide_legend(override.aes = list(size = 4)),  #, #
         colour = guide_legend(override.aes = list(size = 4))) +   #  guide_legend(override.aes = list(size = 5, pch=shapes[1]), order = 2)
  plot_theme(baseSize = 20) + #  size="Year",
  labs( x=glue::glue("C{plot_comps[1]}: {round(percentVar[plot_comps[1]]*100,2)}% variance"),
        y=glue::glue("C{plot_comps[2]}: {round(percentVar[plot_comps[2]]*100,2)}% variance"),
        shape="Year")
# Save plot to a pdf file
ggsave(filedate(glue::glue("Aus_samples_PCA_state_year_PC{paste(plot_comps, collapse='-')}"),
                ext = ".pdf", outdir = "output/plots",
                dateformat = FALSE), width = 10, height=8)
# when subsetting the object, don't forget to subset also the metadata:
# glsub2 <- gl[index.ind, index.loc]
# glsub2@other$ind.metrics <- gl@other$ind.metrics[index.ind,] #not necessary
# glsub2@other$loc.metrics <- gl@other$loc.metrics[index.loc,] #necessary

#### Population Analysis ####
##### Assume one big population #####

# strata_factors <- metadata %>% filter(indNames %in% indNames(glsub)) %>%
#     select(indNames, Location, State, Region, Year, Host, Pathotype, Path_rating, Haplotype, Taxa) %>%
#     mutate(State=factor(State, levels = c("QLD", "NSW", "VIC", "SA", "WA")),
#     Year=factor(Year, levels=as.character(min(Year, na.rm = TRUE):max(Year, na.rm = TRUE))),
#     Host=fct_relevel(fct_infreq(fct_lump_min(Host, min=6)), "Other", after = Inf),
#     Haplotype=fct_relevel(fct_infreq(fct_lump_min(Haplotype, min=2)), "Other", after = Inf),
#     Pathotype=factor(Pathotype, levels=paste0("Group", 0:5))) %>%
#     mutate_at(vars(Year, Host, Haplotype, Pathotype, Region), ~fct_explicit_na(., na_level = "Unknown"))
# glsub@strata <- metadata %>% filter(indNames %in% indNames(glsub)) %>% select(Location, State, Region, Year, Host, Pathotype, Path_rating, Haplotype)
gl_Aus <- glsub
setPop(gl_Aus) <- ~State

#### poppr diversity analysis #####
# convert to genclone object
set.seed(120)
Arab_gclone <- as.genclone(gl2gi(gl_Aus), strata=strata(gl_Aus))
ploidy(Arab_gclone) <- 1
Arab_gclone@strata <- gl_Aus@strata #strata_factors %>% filter(id %in% indNames(Arab_gclone)) # strata(glsub)
# pdf("./output/plots/DArT_genotype_accumulation_curve.pdf", width = 8, height = 6) 
# genotype_curve(Arab_gclone, sample = 1000, quiet = TRUE)
# dev.off()
# assess distance and clustering algorithm
gclone_filtered <- filter_stats(Arab_gclone, plot = TRUE)
thresh_predict <- gclone_filtered %>% map(~cutoff_predictor(.$THRESHOLDS))
alg_choice <- "average" #farthest, nearest
# Contract MLGs based on distance
Arab_dist <- bitwise.dist(Arab_gclone, mat=TRUE, euclidean = TRUE, percent = FALSE)
# xdist <- dist(genind2df(Arab_gclone, usepop = FALSE))
mlg.filter(Arab_gclone, distance = bitwise.dist, percent = FALSE, algorithm = alg_choice) <- thresh_predict[[alg_choice]]

# plot shared MLGs
# splitStrata(monpop) <- ~Tree/Year/Symptom
Arab_tab <- mlg.table(Arab_gclone, strata = ~State, plot = FALSE) %>% as.data.frame() %>% rownames_to_column("State") %>% 
    pivot_longer(-1, names_to = "MLG", values_to = "Count") %>% filter(Count>0) %>% arrange(desc(Count))
# select the top MLGs
shared_MLGs <- Arab_tab %>% count(MLG) %>% arrange(desc(n)) %>% filter(n>1)
MLGs_ordered <- Arab_tab %>% filter(MLG %in% shared_MLGs$MLG) %>% group_by(MLG) %>% summarise(Total=sum(Count)) %>% 
  arrange(desc(Total))
# plot MLGs shared by states
plot_mlgs <- Arab_tab %>% filter(MLG %in% MLGs_ordered$MLG) %>% 
  mutate(State=factor(State, levels=levels(Arab_gclone@strata$State)), MLG=factor(MLG, levels = MLGs_ordered$MLG))
ggplot(plot_mlgs, aes(x=MLG, y=Count, fill=State)) + 
  geom_bar(stat = "identity", width=0.75) + coord_flip() +
  scale_fill_manual(values = my_colours$State)
# save plot
ggsave("output/plots/A_rabiei_poppr_shared_MLG.pdf", width = 7, height=5)
MLG_isolate_map <- Arab_gclone@strata %>% bind_cols(Arab_gclone@mlg@mlg) %>% select(id, State, MLG=contracted)
MLG_isolate_map %>% group_by(MLG) %>% count(State) %>% arrange(desc(n))
Arab_gclone@mlg@mlg %>% count(contracted) %>% arrange(desc(n)) %>% filter(n>1)

# Generate population genetics metrices
region_pop_gen <- poppr(Arab_gclone) %>% mutate(lambda_corr=N/(N - 1) * lambda, CF=MLG/N)  %>% 
  write_xlsx(., "./output/A_rabiei_DArT_poppr.xlsx", "years_by_region", asTable = FALSE, overwritesheet = TRUE)
pdf("./output/plots/A_rabiei_DArT_poppr_rare_state.pdf", width = 8, height = 6) 
region_rare <- diversity_ci(Arab_gclone, n=100, rarefy = TRUE, raw=FALSE) %>% 
  write_xlsx(., "./output/A_rabiei_DArT_poppr.xlsx", "by_region_rare", asTable = FALSE, append=TRUE, overwritesheet = TRUE)
dev.off()

Arab_by_years <- Arab_gclone

Arab_by_years@pop <- Arab_by_years@strata$Year
Arab_by_years <- Arab_by_years[!is.na(pop(Arab_by_years))]

# Generate population genetics metrices
year_pop_gen <- poppr(Arab_by_years) %>% mutate(lambda_corr=N/(N - 1) * lambda, CF=MLG/N) %>%  
  write_xlsx(., "./output/A_rabiei_DArT_poppr.xlsx", "by_year_all_regions", asTable = FALSE, append = TRUE, 
             overwritesheet = TRUE)
pdf("./output/plots/A_rabiei_DArT_poppr_rare_year.pdf", width = 8, height = 6) 
year_pop_gen_rare <- diversity_ci(Arab_by_years, n=100, rarefy = TRUE, raw=FALSE) %>% 
  write_xlsx(., "./output/A_rabiei_DArT_poppr.xlsx", "by_year_rare", asTable = FALSE, append=TRUE, overwritesheet = TRUE)
dev.off()

#### plot haplotype Networks ####
# By state
# setPop(Arab_gclone) <- ~State
Arab_gclone@pop <- Arab_gclone$strata$State
# pe_data <- Arab_gclone # %>%  missingno("genotype", cutoff=0.2)
# my_colours$State
# Calculate the MSN
# pc_colors <- nPop(Arab_gclone) %>% 
#   RColorBrewer::brewer.pal("Set1") %>% 
#   setNames(levels(Arab_gclone@strata$State))
  # setNames(c("QLD", "NSW", "VIC", "SA", "WA"))
# Contract MLGs based on distance
# Arab_gclone@strata <- strata_factors
# setPop(Arab_gclone) <- ~Region
set.seed(120)
mdist <- bitwise.dist(Arab_gclone, mat=TRUE, euclidean = FALSE, percent = FALSE)
pdf(filedate("Arab_regions_msn_all_hosts", ".pdf", here("output/plots")), width = 10, height = 7)
pe_data.msn <- poppr.msn(Arab_gclone, mdist, showplot = TRUE, threshold = thresh_predict[[alg_choice]],
                         palette = my_colours$State, vertex.label = NA, margin=rep(-0.1,4), wscale = FALSE)
dev.off()
# Visualize the network
set.seed(120)
pdf(filedate("Arab_regions_msn_all_hosts2", ".pdf", here("output/plots")), width = 10, height = 7)
plot_poppr_msn(Arab_gclone, pe_data.msn,  palette = my_colours$State,  nodescale = 15, inds = "none",
               margin=rep(0,4))#, nodebase = 1.25) #  
dev.off()

# By Pathotype
Arab_gclone@pop <- Arab_gclone$strata$Patho.Group
# Calculate the MSN
# brewer.pal.info
# pc_colors <-  c("grey", brewer.pal(9, "RdYlGn")[round(seq(1,8, length.out = nPop(Arab_gclone)-1 ) ,0)]) %>% rev(.) %>% 
#   setNames(popNames(Arab_gclone))
# Contract MLGs based on distance
set.seed(120)
pdf(filedate(sprintf("Arab_regions_msn_patho"), ".pdf", "output/plots"), width = 10, height = 7)
pe_data.msn <- poppr.msn(Arab_gclone, mdist, showplot = TRUE, threshold = thresh_predict[[alg_choice]],
                         palette = my_colours$Patho.Group, vertex.label = NA, margin=rep(-0.1,4), wscale = FALSE)
dev.off()

save.image(filedate("A_rabiei_DArT", ".RData"))

# Compare MLG, cluster and WGS analysis ####

mlg_clusters <- cbind(clust_annotation ,MLG=mlg.vector(Arab_gclone), Ind=indNames(Arab_gclone)) %>% 
  arrange(MLG) %>% 
  write_xlsx(., "./output/A_rabiei_DArT_poppr.xlsx", "MLG_Cluster_table", asTable = FALSE, 
             overwritesheet = TRUE)
clusters_in_MLG <- mlg_clusters %>% group_by(MLG) %>% 
  summarise(cluster_num=length(unique(Cluster_fact ))) %>%
  arrange(desc(cluster_num))
mismatch_rate <- nrow(clusters_in_MLG %>% filter(cluster_num>1))/nrow(clusters_in_MLG)

wgs_inds <- readxl::read_excel(here("Sophie_A_rabiei_WGS/A_rabiei_WGS_analysis/sample_info/A_rabiei_isolate_list_for_wgs.xlsx")) %>% select(Isolate) %>% left_join(., mlg_clusters,  by=c("Isolate"="Ind") ) %>% 
  arrange(Cluster)
print(wgs_inds, n=Inf)
missing_inds <- wgs_inds %>% filter(!Isolate %in% mlg_clusters$Ind) %>% .$Isolate 
mlg_clusters$Ind %>% .[grepl("TR64", ., ignore.case = TRUE, perl = TRUE)]
mlg_clusters %>% arrange(Ind)

save.image(filedate("A_rabiei_DArT", ".RData"))

#### Haplotype Analysis ####
# make haplotypes
# dart_haplos <- make_haplotypes(gl3)
# # gl3@ind.names <- ind_names$fixed_names
# ploidy(gl3) <- 1

#### Pathogenicity Association ####
# column_to_rownames("Taxa") %>% mutate(Patho_level=if_else(Path_rating>3, "High", "Low"),
# Patho_fact=factor(Patho_level, levels = c("Low", "High")) )
# SNP.info.pos <- data.frame(snp=locNames(glsub), #,
#                             chr=glsub@other$loc.metrics$Chrom_Ascochyta_rabiei_vX,
#                             pos=glsub@other$loc.metrics$ChromPos_Ascochyta_rabiei_vX+glsub@other$loc.metrics$SnpPosition,
#                             snp_safe_name=gsub("[-\\/]", ".", paste0("X", locNames(glsub))),
#                             stringsAsFactors = FALSE) %>%
# filter(snp_safe_name %in% colnames(export_snps))
# SNP.info.pos <- data.frame(snp=gsub("[-\\/]", ".", paste0("X", locNames(glsub))), #,
#                            chr=glsub@other$loc.metrics$Chrom_Ascochyta_rabiei_vX,
#                            pos=glsub@other$loc.metrics$ChromPos_Ascochyta_rabiei_vX+glsub@other$loc.metrics$SnpPosition,
#                            snp_name=,
#                            stringsAsFactors = FALSE) %>%
#   filter(snp_safe_name %in% colnames(export_snps))
# myData <- gl2sa(glsub) %>% mutate(Isolate=strata_factors$Taxa, Pathogenicity=strata_factors$Path_rating)
myData <- setupSNP(data=export_snps,colSNPs=1:nrow(SNP.info.pos),sep="/",
                    info=SNP.info.pos, sort = TRUE)

export_snps <- gl2sa(glsub) %>%  bind_cols(strata_factors %>% select(Taxa, Path_rating)) %>%
      column_to_rownames("Taxa") %>% mutate(Patho_level=if_else(Path_rating>3, "High", "Low"),
      Patho_fact=factor(Patho_level, levels = c("Low", "High")) )
SNP.info.pos <- data.frame(snp=gsub("[-\\/]", ".", paste0("X", locNames(glsub))), #,
      chr=glsub@other$loc.metrics$Chrom_Ascochyta_rabiei_vX,
      pos=glsub@other$loc.metrics$ChromPos_Ascochyta_rabiei_vX+glsub@other$loc.metrics$SnpPosition,
      snp_name=locNames(glsub),
      stringsAsFactors = FALSE) %>%
      filter(snp %in% colnames(export_snps))
nrow(SNP.info.pos)
SNP.info.pos
# myData <- gl2sa(glsub) %>% mutate(Isolate=strata_factors$Taxa, Pathogenicity=strata_factors$Path_rating)
myData <- setupSNP(data=export_snps,colSNPs=1:nrow(SNP.info.pos),sep="/",
                    info=SNP.info.pos, sort = TRUE)
myData[1:10,1:10]
# Check missingness
# plotMissing(myData)
# Perform association test to Pathogenicity level
res <- WGassociation(Path_rating, data=myData, model="all") # try different models
# res <- scanWGassociation(Pathogenicity, data=myData, model="all")
# Check summary
p_val <- 1e-3
min_snps <- 50
FDR <- 0.05
# Calculate qvalue
res_selected <-res[!is.na(res$`log-additive`)]
res_df <- cbind(attr(res_selected, "gen.info"), "pvalue"=res_selected$`log-additive`) %>%
                bind_cols(., GenABEL::qvaluebh95(.$pvalue, FDR)) %>%
write_xlsx("output/A_rabiei_pathogenicity_snpassoc_results.xlsx", asTable=FALSE, sheet = "snpassoc_results")
# check if it found ny significant SNPs
res_df %>% filter(significant!=FALSE)



