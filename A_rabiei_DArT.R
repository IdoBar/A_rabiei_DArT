# load custom functions from github
devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "util.R")

# install dependencies
# BiocManager::install(c("SNPRelate", "qvalue"))
# install.packages(c("backports", "colorspace", "adegenet"))
# load packages we're going to use
BiocManager::install("qvalue", update = FALSE)
pacman::p_load(  tidyverse, forcats, RColorBrewer, paletteer, ggrepel, glue, dartR, poppr) # quadprog, raster, mvtnorm,

#### Read DArT Data ####
pheno_file <- recent_file("data", "A_rabiei_pathotypes.+.xlsx")
isolate_table <- readxl::read_excel(path = pheno_file , sheet = "pathotyping") %>% select(-one_of(c("Location"))) %>% mutate(Taxa=Isolate) %>%
  mutate(Host=sub("^Genesis09[0-9]+$", "Genesis090", 
                  sub(" ", "", sub("PBA ", "", Host))))
taxa_table <- readxl::read_excel("data/Taxa_DArTseq_Plates1-2_5-05-19.xlsx", na = c("N/A", "")) %>%
  mutate(Host=sub("^Genesis09[0-9]+$", "Genesis090", 
                  sub(" ", "", sub("PBA ", "", Host)))) %>% 
  mutate_at(.vars=c("Year"), as.numeric) %>% select(Taxa, Location,  State, Year, Host, Haplotype) %>% 
  full_join(isolate_table) %>% mutate(State=sub("[^A-Z]", "", toupper(State))) # Host=sub("CICA.+", "CICA1152", Host), 
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

metadata <- tibble(indNames=indNames(gl3)) %>% 
  mutate(Taxa=if_else(is.na(replace_names[indNames]), indNames, replace_names[indNames])) %>% 
  # left_join(isolate_table, by = c("Taxa"="Isolate")) %>% 
  left_join(taxa_table) 
  

metadata %>% 
  select(id=indNames, Location, State, Year, Host, Pathotype, Path_rating, Haplotype)  %>% 
  write_csv("data/DArT_metadata.csv")
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

#### Exploratory Data Analysis ####
gl.report.monomorphs(Arab_gl)  # all loci are polymorphic
# check missing rates
gl.report.callrate(Arab_gl, method = "loc")
gl.report.callrate(Arab_gl, method = "ind", v = 3) # the 2 individuals from Spain are the most remote
# check reproducibility (RepAVG)
gl.report.repavg(Arab_gl)
# check multiple loci per tag and reproducibility
gl.report.secondaries(Arab_gl)
# apply filters
missing_thres <- 0.2
glsub <- Arab_gl %>% gl.filter.callrate(., method = "loc", threshold = (1-missing_thres), v = 3) %>% 
  gl.filter.callrate(., method = "ind", threshold = (1-missing_thres), recalc = TRUE, v = 3) %>% 
  gl.filter.repavg(., threshold = 0.95, v = 3) %>% 
  gl.filter.secondaries(., method = "best", v = 3)

# Convert to allele frequencies (manage missing data)
X <- adegenet::tab(glsub, freq = TRUE, NA.method = "mean")
pca1 <- ade4::dudi.pca(X, scale = FALSE, scannf = FALSE, nf=10)

pca_data <- pca1$li %>% rownames_to_column("indNames")  %>%
  left_join(metadata) %>% mutate(Host=fct_relevel(fct_infreq(fct_lump_min(Host, min=6)), "Other", after = Inf),
                                 State=factor(State, levels = c("QLD", "NSW", "VIC", "SA", "WA", "SPAIN")))
# Calculate variance for each component (columns starting with 'C')
pca_var <- pca1$li %>% summarise_at(vars(starts_with("Axis")), var)
# Calculate the percentage of each component of the total variance
percentVar <- pca_var/sum(pca_var)
# Define colour palette, but get rid of the awful yellow - number 6
pal <- brewer.pal(9, "Set1")
pal2 <- brewer.pal(8, "Dark2")
text_cols <- adjustcolor( c("grey15", "dodgerblue3"), alpha.f = 0.8)
# seq_face <- setNames(c("bold", "plain"),unique(pca_data$Sequenced))
shapes <- c(21:25, 3)
# Isolate to annotate
outliers <- pca_data %>% filter(State=="Spain") 

#### PCA plot 1-2 geographic ####
# shape - Year, Seq.Provider, Host (State)
# Colors - Haplotype (qualitative), State (qualitative), Pathotype (diverging - RdYlGn)
# Size - Pathogenicity
# Create the plot for C1 and C2
ggplot(pca_data, aes(x=Axis1, y=Axis2, fill=State)) + 
  geom_point(alpha = 0.85, size=5, shape=21) +
  scale_fill_paletteer_d(ggsci, category10_d3 ) +
  plot_theme(baseSize = 20) + #  size="Year",
  labs( x=glue::glue("C1: {round(percentVar[1]*100,2)}% variance"),
        y=glue::glue("C2: {round(percentVar[2]*100,2)}% variance"))
# Save plot to a pdf file
ggsave(filedate(glue::glue("EDA_all_samples_PCA_state_PC1-2"),
                ext = ".pdf", outdir = "output/plots",
                dateformat = FALSE), width = 10, height=8)


#### Australian Pops Analysis ####
# filter on missing data (will remove Spanish isolates)
# apply filters
missing_thres <- 0.1
glsub <- Arab_gl %>% gl.filter.callrate(., method = "loc", threshold = (1-missing_thres), v = 3) %>% 
  gl.filter.callrate(., method = "ind", threshold = (1-missing_thres), recalc = TRUE, v = 3) %>% 
  gl.filter.repavg(., threshold = 0.95, v = 3) %>% 
  gl.filter.secondaries(., method = "best", v = 3)

# Convert to allele frequencies (manage missing data)
X <- adegenet::tab(glsub, freq = TRUE, NA.method = "mean")
pca1 <- ade4::dudi.pca(X, scale = FALSE, scannf = FALSE, nf=10)

pca_data <- pca1$li %>% rownames_to_column("indNames")  %>%
  left_join(metadata) %>% mutate(Host=fct_relevel(fct_infreq(fct_lump_min(Host, min=6)), "Other", after = Inf),
                                 State=factor(State, levels = c("QLD", "NSW", "VIC", "SA", "WA")))
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
outliers <- pca_data %>% filter(State=="Spain") 
# outliers <- c("16CUR017", "TR9543", "TR9529", "15CUR003") # PacBio-sequenced
# outliers <- pca_data$Isolate
#### PCA plots biological ####
# shape - Year, Seq.Provider, Host (State)
# Colors - Haplotype (qualitative), State (qualitative), Pathotype (diverging - RdYlGn)
# Size - Pathogenicity
# Create the plot for C1 and C2
plot_comps <- 1:2
plot_axes <- paste0("Axis", plot_comps)
ggplot(pca_data, aes(x=!!sym(plot_axes[1]), y=!!sym(plot_axes[2]), shape=State, fill=Pathotype)) + 
  geom_point(alpha = 0.85, size=5) +
  scale_fill_brewer(palette = "RdYlGn", direction = -1) +
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
                ext = ".pdf", outdir = "output/plots",
                dateformat = FALSE), width = 10, height=8)
# Create the plot for C3 and C4
plot_comps <- 3:4
plot_axes <- paste0("Axis", plot_comps)
ggplot(pca_data, aes(x=!!sym(plot_axes[1]), y=!!sym(plot_axes[2]), shape=State, fill=Pathotype)) + 
  geom_point(alpha = 0.85, size=5) +
  scale_fill_brewer(palette = "RdYlGn", direction = -1) + 
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
                ext = ".pdf", outdir = "output/plots",
                dateformat = FALSE), width = 10, height=8)


#### Highlight high pathogenicity isolates  ####
# Create the plot for C1 and C2
plot_comps <- 3:4
plot_axes <- paste0("Axis", plot_comps)
path_cols <- c(rep(NA, 4),"#FDAE61" ,"#D73027", NA)

ggplot(pca_data, aes(x=!!sym(plot_axes[1]), y=!!sym(plot_axes[2]), shape=State, fill=Pathotype)) + 
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
                ext = ".pdf", outdir = "output/plots",
                dateformat = FALSE), width = 10, height=8)

#### PCA plots spatio-temporal ####
# shape - Year, Seq.Provider, Host (State)
# Colors - Haplotype (qualitative), State (qualitative), Pathotype (diverging - RdYlGn)
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



#### Haplotype Analysis ####
# make haplotypes
dart_haplos <- make_haplotypes(gl3)
# gl3@ind.names <- ind_names$fixed_names
ploidy(gl3) <- 1


as.matrix(gl3)[,1:10]
