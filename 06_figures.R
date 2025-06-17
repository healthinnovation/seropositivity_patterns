library(tidyverse)
library(ComplexUpset)
library(ggnewscale)

ffi <- read_csv("./data/mid/ffi_hh_ind_seropos.csv")


# Upset Plot --------------------------------------------------------------

# Data Preparation
df_bin <- ffi %>% 
  mutate(
    schistosomiasis = if_else(
      SEA_pos == 1 | X229.E.NP_pos == 1,
      1,
      0
    ),
    covid = if_else(
      MERSP.NP_pos == 1 | SARS.NP.WT_pos == 1 | SARS.RBP.WT_pos == 1,
      1,
      0
    ),
    chik = if_else(
      Chik.E1_pos == 1, 
      1, 
      0
    ),
    zika = if_else(
      Zika.NS1_pos == 1,
      1,
      0
    ),
    malaria = if_else(
      pv_exposure == "Positive" | pf_exposure == "Positive",
      1,
      0
    )
  )%>% 
  select(ffi_is_code, schistosomiasis, covid, chik, zika,
         malaria) %>% 
  drop_na()


# Plot

upset_plot <- ComplexUpset::upset(
  df_bin,
  intersect = c(
    "schistosomiasis",
    "covid",
    "chik",
    "zika",
    "malaria"
  ),
  name = "Exposure Profile",
  base_annotations = list(
    'Count' = intersection_size(
      mapping = aes(fill = after_stat(y))
    ) +
      scale_fill_gradient(
        low = "#e6f2f2",
        high = "#294d61",
        guide = "none") +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
      )
  ),
  queries = list(
    upset_query(set = "covid", fill = "#036666"),
    upset_query(set = "zika", fill = "#14746f"),
    upset_query(set = "chik", fill = "#469d89"),
    upset_query(set = "schistosomiasis", fill = "#67b99a"),
    upset_query(set = "malaria", fill = "#78c6a3")
  ),
  set_sizes = upset_set_size() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    ) +
    ylab("Set")
)

upset_plot

# Spatial Analysis --------------------------------------------------------

df_spat <- ffi %>% 
  mutate(
    schistosomiasis = if_else(
      SEA_pos == 1 | X229.E.NP_pos == 1,
      1,
      0
    ),
    covid = if_else(
      MERSP.NP_pos == 1 | SARS.NP.WT_pos == 1 | SARS.RBP.WT_pos == 1,
      1,
      0
    ),
    chik = if_else(
      Chik.E1_pos == 1, 
      1, 
      0
    ),
    zika = if_else(
      Zika.NS1_pos == 1,
      1,
      0
    ),
    malaria = if_else(
      pv_exposure == "Positive" | pf_exposure == "Positive",
      1,
      0
    )
  ) %>% 
  select()







