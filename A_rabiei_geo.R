# load custom functions from github
devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "util.R")

# install dependencies
# gganimate requires gifski to render gifs, which first requires rust, see installation instructions: https://github.com/r-rust/hellorust#installation
# or download from this link: https://gif.ski/ and put gifski.exe in a folder in your Path
# devtools::install_github("r-rust/gifski")
# load packages we're going to use
# BiocManager::install("qvalue", update = FALSE)
pacman::p_load(  tidyverse, forcats, paletteer, ggrepel, glue, ggmap, gganimate, ggimage, gifski) # quadprog, raster, mvtnorm, dartR

# register Google API key
register_google(key = "AIzaSyCdLpTphIZnduo2g7XDCJjB6W5Z0qtBau0", write = TRUE)

#### Read DArT Data ####
pheno_file <- recent_file("data", "A_rabiei_pathotypes.+.xlsx")
isolate_table <- readxl::read_excel(path = pheno_file , sheet = "pathotyping") %>% select(-one_of(c("Location"))) %>% mutate(Taxa=Isolate) %>%
  mutate(Host=sub("^Genesis09[0-9]+$", "Genesis090", 
                  sub(" ", "", sub("PBA ", "", Host))))
taxa_table <- readxl::read_excel("data/Taxa_DArTseq_Plates1-2_5-05-19.xlsx", na = c("N/A", "")) %>%
  mutate(Host=sub("^Genesis09[0-9]+$", "Genesis090", 
                  sub(" ", "", sub("PBA ", "", Host)))) %>% 
  mutate_at(.vars=c("Year"), as.numeric) %>% select(Taxa, Location,  State, Year, Host, Haplotype) %>% 
  write_xlsx(., excel_file = "output/A_rabiei_SSR_locations.xlsx", sheet = "Taxa_locations")
  full_join(isolate_table) %>% mutate(State=sub("[^A-Z]", "", toupper(State))) 

taxa_table %>% # Host=sub("CICA.+", "CICA1152", Host), 
  write_xlsx("data/5_year_complete_SSR_db.xlsx", "2013-2018_pheno", append=TRUE, asTable=FALSE)
ssr_haplos <- readxl::read_excel("data/5_year_complete_SSR_db.xlsx")

taxa_table %>% mutate(Address=paste(Location, State, sep=", ")) %>% group_by(Address, Year) %>% 
  count(Address) %>% write_xlsx(., excel_file = "output/A_rabiei_SSR_locations.xlsx", sheet = "Taxa_locations")

ssr_location_sum <- readxl::read_excel("output/A_rabiei_SSR_locations_curated.xlsx") %>% group_by(Address, Year) %>% 
  summarise(appearances=sum(n)) 

# summarise DArT samples by year and location
dart_sum <- read_csv("data/DArT_metadata.csv") %>% filter(!is.na(Year), State!="SPAIN") %>% 
  group_by(State, Year) %>% 
  count() %>% pivot_wider(names_from = State, values_from = n) %>% 
  ungroup() %>% mutate_all(~replace_na(.x, 0)) %>% arrange(Year) %>% 
  write_xlsx(., "output/dart_samples_sum.xlsx", sheet = "sheet1")

source("src/nominatim_osm.R")
ssr_adds <- purrr::map_dfr(unique(ssr_location_sum$Address), nominatim_osm) 

ssr_locations <- inner_join(ssr_location_sum, ssr_adds, by=c("Address"="address")) %>% filter(Year>2010, lon>100, lat < (-20))
ssr_locations %>% arrange(desc(lat))


# plot map
shapes <- c(21:25)
aus <- c(left= 111, bottom = -40, right = 156.5, top = -13.4)
collection_map <- get_stamenmap(aus, zoom = 5, maptype = "terrain") %>% 
  ggmap(base_layer = ggplot(ssr_locations, 
                            mapping = aes(x = lon, y = lat,  size = appearances, group=factor(Year), shape=factor(Year), fill=factor(Year)))) +
  geom_point(alpha = .65) +
  guides(shape = guide_legend(override.aes = list(size = 4)), size = guide_legend(order = 2)) +
  scale_size_continuous(range = c(3, 8), breaks = c(5,10,20,40)) + # colour = 'blue',
  scale_shape_manual(values = shapes) + 
  scale_fill_paletteer_d(ggsci, category10_d3 ) +
  labs(size = 'Samples', title = 'Samples collected in {round(frame_time,0)}', fill="Year", shape="Year", group="Year") + plot_theme(baseSize = 20) 
  

ggsave("output/plots/Arab_geo_collection.pdf", width = 12, height = 10)
# animate the plot
animate(collection_map +
          transition_time(Year) + enter_fade() + shadow_mark() #+ ease_aes('cubic-in-out') , renderer = magick_renderer()
  , width = 800, height = 600, duration = 20, end_pause = 10)

anim_save("output/sample_collection_regions_by_year_gifski.gif")

#### Plot production ####
# prepare ABARES data
# download from http://data.daff.gov.au/data/warehouse/aucrpd9abcc003/aucrpd9aba_20180911_Z0Srg/AustCropRrt20180911_StateCropData_v1.0.0.xlsx and edit to combine data from all states
states <- set_names(c("QLD", "NSW", "VIC", "SA", "WA"), c("Queensland", "New South Wales", "Victoria", "South Australia", "Western Australia"))
prod_data <- readxl::read_excel("data/Chickpea_production_2003_2019.xlsx", sheet = "Chickpea_prod") %>% 
  tidyr::pivot_longer(starts_with("20"), names_to = "Period", values_to = "value") %>% 
  tidyr::pivot_wider(names_from = Measure, values_from = value) %>% 
  mutate(Period=sub(" (\\w)$", "(\\1)", Period)) %>% 
  mutate(State_abbr=factor(states[State], levels = states), State=factor(State, levels = names(states)))
prod_data %>% group_by(Period) %>% summarise(sum(Production))

# plot stacked bar graph
ggplot(prod_data, aes(x=Period, y=Production, fill=State_abbr, group=State_abbr)) +
  # geom_line(size=1) +
  geom_bar(stat = "identity", width = 0.5) +
  # geom_line(size=1.5) +
  scale_fill_paletteer_d(RColorBrewer, Set1) + scale_y_continuous(expand = expand_scale(add = c(50, 200))) + # ggsci, category10_d3 ; ggthemes, wsj_colors6 
  # coord_cartesian(ylim = c(-25, 2200)) + 
  labs(x="ABARES Reporting Period", y="Production (kt)", fill="State") +
  # coord_flip(expand = FALSE, ylim = c(-25, 2200)) +
  plot_theme() + theme(axis.text.x = element_text(angle = -45, hjust = 0))
  
ggsave("output/plots/Chickpea_production_ABARES_abbr.pdf", height = 8, width = 12)

# plot line graph
ggplot(prod_data, aes(x=Period, y=Production, colour=State, fill=State, group=State, shape=State)) +
  geom_line(size=1) + geom_point(size=4) + scale_shape_manual(values=21:25) +
  scale_colour_paletteer_d(RColorBrewer, Set1) +
  scale_fill_paletteer_d(RColorBrewer, Set1) + scale_y_continuous(expand = expand_scale(add = c(50, 200))) + # ggsci, category10_d3 ; ggthemes, wsj_colors6 
  # coord_cartesian(ylim = c(-25, 2200)) + 
  labs(x="ABARES Reporting Period", y="Production (kt)") +
  # coord_flip(expand = FALSE, ylim = c(-25, 2200)) +
  plot_theme() + theme(axis.text.x = element_text(angle = -45, hjust = 0))

ggsave("output/plots/Chickpea_production_ABARES_line.pdf", height = 8, width = 12)

price_data <- readxl::read_excel("data/Chickpea_production_2003_2019.xlsx", sheet = "Chickpea_prices")
# plot line graph
ggplot(price_data, aes(x=Quarter, y=Price, group=1)) +
  geom_line(size=4, colour="grey30")  + #+ geom_point(size=4, colour="grey10", pch=12) 
  # coord_cartesian(ylim = c(-25, 2200)) + 
  labs(x="ABARES Reporting Period", y="Price (A$/t)") +
  # coord_flip(expand = FALSE, ylim = c(-25, 2200)) +
  plot_theme() + theme(axis.text.x = element_text(angle = -45, hjust = 0), panel.grid = element_blank(), panel.background = element_blank())

ggsave("output/plots/Chickpea_prices_ABARES_line.pdf", height = 8, width = 12)
