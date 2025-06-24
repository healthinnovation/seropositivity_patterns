library(tidyverse)
library(ComplexUpset)
library(sf)
library(sfdep)
library(ggspatial)
library(patchwork)

ffi <- read_csv("./data/mid/ffi_hh_ind_seropos.csv")


# Upset Plot --------------------------------------------------------------

# Data Preparation
df_bin <- ffi %>% 
  mutate(
    Schistosomiasis = if_else(
      SEA_pos == 1 | X229.E.NP_pos == 1,
      1,
      0
    ),
    Covid = if_else(
      MERSP.NP_pos == 1 | SARS.NP.WT_pos == 1 | SARS.RBP.WT_pos == 1,
      1,
      0
    ),
    Chik = if_else(
      Chik.E1_pos == 1, 
      1, 
      0
    ),
    Zika = if_else(
      Zika.NS1_pos == 1,
      1,
      0
    ),
    Malaria = if_else(
      pv_exposure == "Positive" | pf_exposure == "Positive",
      1,
      0
    )
  )%>% 
  select(ffi_is_code, Schistosomiasis, Covid, Chik, Zika,
         Malaria) %>% 
  drop_na()


# Plot

upset_plot <- ComplexUpset::upset(
  df_bin,
  intersect = c(
    "Schistosomiasis",
    "Covid",
    "Chik",
    "Zika",
    "Malaria"
  ),
  name = "Exposure Profile",
  base_annotations = list(
    'Count' = intersection_size(
      mapping = aes(fill = after_stat(y))
    ) +
      scale_fill_gradient(
        low = "#a1cca5",
        high = "#415d43",
        guide = "none") +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
      )
  ),
  queries = list(
    upset_query(set = "Covid", fill = "#344e41"),
    upset_query(set = "Zika", fill = "#3a5a40"),
    upset_query(set = "Chik", fill = "#588157"),
    upset_query(set = "Schistosomiasis", fill = "#a3b18a"),
    upset_query(set = "Malaria", fill = "#dad7cd")
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

ggsave("output/upset_plot.png",
       upset_plot,
       dpi = 1200,
       width = 13,
       height = 8,
       device = grDevices::png)


# Spatial Analysis --------------------------------------------------------

# Data Preparation 
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
  select(ffi_h_code, ffi_is_code,
         schistosomiasis, covid,
         chik, zika, malaria,
         ffi_gps_lat, ffi_gps_long) %>% 
  group_by(ffi_h_code) %>%
  summarise(
    prev_schisto  = sum(schistosomiasis, na.rm = TRUE) / n(),   
    prev_covid    = sum(covid,            na.rm = TRUE) / n(),
    prev_chik     = sum(chik,             na.rm = TRUE) / n(),
    prev_zika     = sum(zika,             na.rm = TRUE) / n(),
    prev_malaria  = sum(malaria,          na.rm = TRUE) / n(),
    lat  = first(ffi_gps_lat),
    long = first(ffi_gps_long)
  ) %>% 
  st_as_sf(
    coords = c("long", "lat"),
    crs = 4326
  )


# Neighbours and spatial weights 
df_n <- df_spat %>% 
  mutate(
    nb = st_knn(geometry, k = 5),
    wt = st_weights(nb)
  )


# Spatial: Schistosomiasis
local_m_schisto <- df_n %>% 
  mutate(
    lm_schisto = local_moran(prev_schisto, nb, wt)
  ) %>% 
  unnest(lm_schisto)

# Mapa
data("Peru", package = "innovar")

crop_coordinates <- Peru %>% 
  filter(dep == "LORETO",
         distr %in% c("INDIANA",
                      "BELEN"))


map_schisto <- local_m_schisto %>% 
  mutate(mean_cat = if_else(p_folded_sim <= 0.1, as.character(mean),
                            "Not Significant"),
         mean_cat = factor(mean_cat,
                           levels = c(
                             "High-High",
                             "Low-High",
                             "High-Low",
                             "Low-Low",
                             "Not Significant"
                           ))) %>% 
  arrange(mean_cat != "Not Significant") %>% 
  ggplot(aes(geometry = geometry, colour = mean_cat)) +
  annotation_map_tile(type = "cartolight", zoom = 12) +
  geom_sf(data = crop_coordinates,
          fill = NA,
          color = "black",
          linetype = "dashed",
          linewidth = 0.6) +
  geom_sf(size = 3, alpha = 0.9) +
  scale_colour_manual(
    values = c(
      "Not Significant" = "grey40",
      "Low-Low" = "#03595c",
      "High-Low" = "#87cdb5",
      "Low-High" = "#efc38f",
      "High-High" = "#e2634b"
    )
  ) +
  annotation_scale() +
  annotation_north_arrow(location = "tl", style = north_arrow_nautical) +
  theme_bw() +
  labs(
    colour = "Cluster"
  ) 


# Spatial: Malaria
local_m_malaria <- df_n %>% 
  mutate(
    lm_malaria = local_moran(prev_malaria, nb, wt)
  ) %>% 
  unnest(lm_malaria)

# Mapa

map_malaria <- local_m_malaria %>% 
  mutate(mean_cat = if_else(p_folded_sim <= 0.1, as.character(mean),
                            "Not Significant"),
         mean_cat = factor(mean_cat,
                           levels = c(
                             "High-High",
                             "Low-High",
                             "High-Low",
                             "Low-Low",
                             "Not Significant"
                           ))) %>% 
  arrange(mean_cat != "Not Significant") %>% 
  ggplot(aes(geometry = geometry, colour = mean_cat)) +
  annotation_map_tile(type = "cartolight", zoom = 12) +
  geom_sf(data = crop_coordinates,
          fill = NA,
          color = "black",
          linetype = "dashed",
          linewidth = 0.6) +
  geom_sf(size = 3, alpha = 0.9) +
  scale_colour_manual(
    values = c(
      "Not Significant" = "grey40",
      "Low-Low" = "#03595c",
      "High-Low" = "#87cdb5",
      "Low-High" = "#efc38f",
      "High-High" = "#e2634b"
    )
  ) +
  annotation_scale() +
  annotation_north_arrow(location = "tl", style = north_arrow_nautical) +
  theme_bw() +
  labs(
    colour = "Cluster"
  ) +
  theme(
    legend.position = "none"
  )


mal_schisto <- (map_malaria + map_schisto) +
  plot_annotation(
    tag_levels = c("a", "b"),
    tag_suffix = ")",
    theme = theme(
      plot.tag = element_text(
        size = 16, face = "bold"
      ),
      plot.tag.position = c(0.02, 0.98)
    )
  )

ggsave("output/spatial_malaria_schisto.png",
       mal_schisto,
       dpi = 600,
       width = 11,
       height = 6,
       device = grDevices::png)




# Spatial: Zika
local_m_zika <- df_n %>% 
  mutate(
    lm_zika = local_moran(prev_zika, nb, wt)
  ) %>% 
  unnest(lm_zika)

# Mapa

map_zika <- local_m_zika %>% 
  mutate(mean_cat = if_else(p_folded_sim <= 0.1, as.character(mean),
                            "Not Significant"),
         mean_cat = factor(mean_cat,
                           levels = c(
                             "High-High",
                             "Low-High",
                             "High-Low",
                             "Low-Low",
                             "Not Significant"
                           ))) %>% 
  arrange(mean_cat != "Not Significant") %>%
  ggplot(aes(geometry = geometry, colour = mean_cat)) +
  annotation_map_tile(type = "cartolight", zoom = 12) +
  geom_sf(data = crop_coordinates,
          fill = NA,
          color = "black",
          linetype = "dashed",
          linewidth = 0.6) +
  geom_sf(size = 3, alpha = 0.9) +
  scale_colour_manual(
    values = c(
      "Not Significant" = "grey40",
      "Low-Low" = "#03595c",
      "High-Low" = "#87cdb5",
      "Low-High" = "#efc38f",
      "High-High" = "#e2634b"
    )
  ) +
  annotation_scale() +
  annotation_north_arrow(location = "tl", style = north_arrow_nautical) +
  theme_bw() +
  labs(
    colour = "Cluster"
  )



# Spatial: Chik
local_m_chik <- df_n %>% 
  mutate(
    lm_chik = local_moran(prev_chik, nb, wt)
  ) %>% 
  unnest(lm_chik)

# Mapa

map_chik <- local_m_chik %>% 
  mutate(mean_cat = if_else(p_folded_sim <= 0.1, as.character(mean),
                            "Not Significant"),
         mean_cat = factor(mean_cat,
                           levels = c(
                             "High-High",
                             "Low-High",
                             "High-Low",
                             "Low-Low",
                             "Not Significant"
                           ))) %>% 
  arrange(mean_cat != "Not Significant") %>%
  ggplot(aes(geometry = geometry, colour = mean)) +
  annotation_map_tile(type = "cartolight", zoom = 12) +
  geom_sf(data = crop_coordinates,
          fill = NA,
          color = "black",
          linetype = "dashed",
          linewidth = 0.6) +
  geom_sf(size = 3, alpha = 0.9) +
  scale_colour_manual(
    values = c(
      "Not Significant" = "grey40",
      "Low-Low" = "#03595c",
      "High-Low" = "#87cdb5",
      "Low-High" = "#efc38f",
      "High-High" = "#e2634b"
    )
  ) +
  annotation_scale() +
  annotation_north_arrow(location = "tl", style = north_arrow_nautical) +
  theme_bw() +
  labs(
    colour = "Cluster"
  ) +
  theme(
    legend.position = "none"
  )

zik_chik <- map_chik + map_zika + 
  plot_annotation(
    tag_levels = "a",
    tag_suffix = ")",
    theme = theme(
      plot.tag = element_text(
        size = 16, face = "bold"
      ),
      plot.tag.position = c(0.02, 0.98)
    )
  ) 

zik_chik

ggsave("output/spatial_zik_chik.png",
       zik_chik,
       dpi = 600,
       width = 11,
       height = 6,
       device = grDevices::png)


# Spatial: Covid
local_m_covid <- df_n %>% 
  mutate(
    lm_covid = local_moran(prev_covid, nb, wt)
  ) %>% 
  unnest(lm_covid)

# Mapa
map_covid <- local_m_covid %>% 
  mutate(mean_cat = if_else(p_folded_sim <= 0.1, as.character(mean),
                            "Not Significant"),
         mean_cat = factor(mean_cat,
                           levels = c(
                             "High-High",
                             "Low-High",
                             "High-Low",
                             "Low-Low",
                             "Not Significant"
                           ))) %>% 
  arrange(mean_cat != "Not Significant") %>%
  ggplot(aes(geometry = geometry, colour = mean_cat)) +
  annotation_map_tile(type = "cartolight", zoom = 12) +
  geom_sf(data = crop_coordinates,
          fill = NA,
          color = "black",
          linetype = "dashed",
          linewidth = 0.6) +
  geom_sf(size = 3, alpha = 0.9) +
  scale_colour_manual(
    values = c(
      "Not Significant" = "grey40",
      "Low-Low" = "#03595c",
      "High-Low" = "#87cdb5",
      "Low-High" = "#efc38f",
      "High-High" = "#e2634b"
    )
  ) +
  annotation_scale() +
  annotation_north_arrow(location = "tl", style = north_arrow_nautical) +
  theme_bw() +
  labs(
    colour = "Cluster"
  )

ggsave("output/spatial_covid.png",
       map_covid,
       dpi = 600,
       width = 8,
       height = 7,
       device = grDevices::png)
