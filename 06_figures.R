library(tidyverse)
library(ComplexUpset)
library(sf)
library(sfdep)
library(ggspatial)
library(patchwork)
library(prettymapr)
library(osmdata)

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
  ) %>% 
  select(ffi_is_code, Schistosomiasis, Chik, Zika, Malaria) %>% 
  drop_na()

# Plot
upset_plot <- ComplexUpset::upset(
  df_bin,
  intersect = c("Schistosomiasis", "Chik", "Zika", "Malaria"),
  name = "Exposure Profile",
  base_annotations = list(
    'Count' = intersection_size(
      mapping = aes(fill = after_stat(y))
    ) +
      scale_fill_gradient(
        low = "#67a9cf",
        high = "#67a9cf",
        guide = "none"
      ) +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
      )
  ),
  queries = list(
    upset_query(set = "Zika", fill = "#67a9cf"),
    upset_query(set = "Chik", fill = "#67a9cf"),
    upset_query(set = "Schistosomiasis", fill = "#67a9cf"),
    upset_query(set = "Malaria", fill = "#67a9cf")
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

# Mostrar el gráfico
upset_plot

# Guardar como imagen
ggsave("output/upset_plot.png",
       upset_plot,
       dpi = 1200,
       width = 13,
       height = 8)


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

local_m_schisto_rot <- local_m_schisto |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')
crop_coordinates_rot <- crop_coordinates |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')

bbox_crop <- st_bbox(c(xmin = -3179322 , xmax = 3199161,
                       ymin = -3696615 , ymax = -3603883),
                     crs = st_crs(local_m_schisto_rot))

local_m_schisto_crop <- st_crop(local_m_schisto_rot, bbox_crop)

map_schisto <- local_m_schisto_crop %>%
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
  annotation_map_tile(type = "osm", zoom = 11) +
  geom_sf(data = crop_coordinates,
          fill = NA,
          color = "black",
          linetype = "dashed",
          linewidth = 0.6) +
  geom_sf(size = 3, alpha = 0.8) +
  scale_colour_manual(
    values = c(
      "Not Significant" = "grey40",
      "Low-Low" = "#2166ac",
      "High-Low" = "#67a9cf",
      "Low-High" = "#ef8a62",
      "High-High" = "#b2182b"
    )
  ) +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "tl", style = north_arrow_nautical, rotation = +40) +
  theme_bw() +
  labs(
    colour = "Cluster"
    ) + theme(
      legend.position = "none"
    ) +
  coord_sf(
    xlim = c(3172322, 3199761),
    ylim = c(-3698615, -3603883)
  )


# Spatial: Malaria
local_m_malaria <- df_n %>% 
  mutate(
    lm_malaria = local_moran(prev_malaria, nb, wt)
  ) %>% 
  unnest(lm_malaria)

local_m_malaria_rot <- local_m_malaria |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')

bbox_crop <- st_bbox(c(xmin = -3179322 , xmax = 3199161,
                       ymin = -3696615 , ymax = -3603883),
                     crs = st_crs(local_m_malaria_rot))

local_m_malaria_crop <- st_crop(local_m_malaria_rot, bbox_crop)

# Mapa

map_malaria <- local_m_malaria_crop %>% 
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
  annotation_map_tile(type = "osm", zoom = 11) +
  geom_sf(data = crop_coordinates,
          fill = NA,
          color = "black",
          linetype = "dashed",
          linewidth = 0.6) +
  geom_sf(size = 3, alpha = 0.8) +
  scale_colour_manual(
    values = c(
      "Not Significant" = "grey40",
      "Low-Low" = "#2166ac",
      "High-Low" = "#67a9cf",
      "Low-High" = "#ef8a62",
      "High-High" = "#b2182b"
    )
  ) +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "tl", style = north_arrow_nautical, rotation = +40) +
  theme_bw() +
  labs(
    colour = "Cluster"
  ) + theme(
    legend.position = "none"
  ) +
  coord_sf(
    xlim = c(3172322, 3199761),
    ylim = c(-3698615, -3603883)
  )


# Spatial: Zika
local_m_zika <- df_n %>% 
  mutate(
    lm_zika = local_moran(prev_zika, nb, wt)
  ) %>% 
  unnest(lm_zika)

local_m_zika_rot <- local_m_zika |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')

bbox_crop <- st_bbox(c(xmin = -3179322 , xmax = 3199161,
                       ymin = -3696615 , ymax = -3603883),
                     crs = st_crs(local_m_zika_rot))

local_m_zika_crop <- st_crop(local_m_zika_rot, bbox_crop)

# Mapa

map_zika <- local_m_zika_crop %>% 
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
  annotation_map_tile(type = "osm", zoom = 11) +
  geom_sf(data = crop_coordinates,
          fill = NA,
          color = "black",
          linetype = "dashed",
          linewidth = 0.6) +
  geom_sf(size = 3, alpha = 0.8) +
  scale_colour_manual(
    values = c(
      "Not Significant" = "grey40",
      "Low-Low" = "#2166ac",
      "High-Low" = "#67a9cf",
      "Low-High" = "#ef8a62",
      "High-High" = "#b2182b"
    )
  ) +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "tl", style = north_arrow_nautical, rotation = +40) +
  theme_bw() +
  labs(
    colour = "Cluster"
  ) +
  coord_sf(
    xlim = c(3172322, 3199761),
    ylim = c(-3698615, -3603883)
  )



# Spatial: Chik
local_m_chik <- df_n %>% 
  mutate(
    lm_chik = local_moran(prev_chik, nb, wt)
  ) %>% 
  unnest(lm_chik)

local_m_chik_rot <- local_m_chik |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')

bbox_crop <- st_bbox(c(xmin = -3179322 , xmax = 3199161,
                       ymin = -3696615 , ymax = -3603883),
                     crs = st_crs(local_m_chik_rot))

local_m_chik_crop <- st_crop(local_m_chik_rot, bbox_crop)

# Mapa

map_chik <- local_m_chik_crop %>% 
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
  annotation_map_tile(type = "osm", zoom = 11) +
  geom_sf(data = crop_coordinates,
          fill = NA,
          color = "black",
          linetype = "dashed",
          linewidth = 0.6) +
  geom_sf(size = 3, alpha = 0.8) +
  scale_colour_manual(
    values = c(
      "Not Significant" = "grey40",
      "Low-Low" = "#2166ac",
      "High-Low" = "#67a9cf",
      "Low-High" = "#ef8a62",
      "High-High" = "#b2182b"
    )
  ) +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "tl", style = north_arrow_nautical, rotation = +40) +
  theme_bw() +
  labs(
    colour = "Cluster"
  ) + theme(
    legend.position = "none"
  ) +
  coord_sf(
    xlim = c(3172322, 3199761),
    ylim = c(-3698615, -3603883)
  )

# Joining maps

map_four_join <- map_schisto + map_malaria + map_chik + map_zika +
  plot_layout(
    ncol = 4
  ) +
  plot_annotation(
    tag_levels = "a",
    tag_suffix = ")",
    theme = theme(
      plot.tag = element_text(
        size=16, face = "bold"
      ),
      plot.tag.position = c(0.02, 0.98)
    )
  )

ggsave("output/spatial_combine_four.png",
       map_four_join,
       dpi = 600,
       width = 14,
       height = 7.5,
       device = grDevices::png)



# Omit this code -------------------------------------------------------------------------------

# Spatial: Covid
local_m_covid <- df_n %>% 
  mutate(
    lm_covid = local_moran(prev_covid, nb, wt)
  ) %>% 
  unnest(lm_covid)

local_m_covid_rot <- local_m_covid |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')

bbox_crop <- st_bbox(c(xmin = -3179322 , xmax = 3199161,
                       ymin = -3696615 , ymax = -3603883),
                     crs = st_crs(local_m_covid_rot))

local_m_covid_crop <- st_crop(local_m_covid_rot, bbox_crop)

# Mapa
map_covid <- local_m_covid_crop %>% 
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
  geom_sf(size = 3, alpha = 0.8) +
  scale_colour_manual(
    values = c(
      "Not Significant" = "grey40",
      "Low-Low" = "#03595c",
      "High-Low" = "#87cdb5",
      "Low-High" = "#efc38f",
      "High-High" = "#e2634b"
    )
  ) +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "tl", style = north_arrow_nautical, rotation = +40) +
  theme_bw() +
  labs(
    colour = "Cluster"
  ) + theme(
    legend.position = "none"
  ) +
  coord_sf(
    xlim = c(3172322, 3199761),
    ylim = c(-3698615, -3603883)
  )
