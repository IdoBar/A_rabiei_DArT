# load custom functions from github
devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "util.R")

# install dependencies
# gganimate requires gifski to render gifs, which first requires rust, see installation instructions: https://github.com/r-rust/hellorust#installation
# or download from this link: https://gif.ski/ and put gifski.exe in a folder in your Path
# devtools::install_github("r-rust/gifski")
# load packages we're going to use
# BiocManager::install("qvalue", update = FALSE)
pacman::p_load(tidyverse, paletteer) # quadprog, raster, mvtnorm, dartR


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
