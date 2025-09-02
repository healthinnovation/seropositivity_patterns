library(tidyverse)
library(ComplexUpset)
library(sf)
library(sfdep)
library(ggspatial)
library(patchwork)
library(prettymapr)
library(osmdata)

ffi <- read_csv("./data/mid/ffi_hh_ind_seropos.csv") %>% 
  mutate(age_oms = case_when(
    age_cat %in% c("[0-10)", "[10-20)") ~ "Children and adolescent",
    age_cat %in% c("[20-30)", "[30-40)") ~ "Young adult",
    age_cat %in% c("[40-50)", "[50-60)") ~ "Middle-aged adult",
    age_cat %in% c("[60-70)", "[70+)")   ~ "Older adult",
    TRUE ~ NA_character_
    )
    )


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
    Pvivax = if_else(
      pv_exposure == "Positive",
      1,
      0
    ),
    Pfalciparum = if_else(
      pf_exposure == "Positive",
      1,
      0
    )
  ) %>% 
  select(ffi_is_code, Schistosomiasis, Chik, Zika, Pvivax,
         Pfalciparum, age_oms) %>% 
  drop_na()

# Plot
upset_plot <- ComplexUpset::upset(
  df_bin,
  intersect = c("Schistosomiasis", "Chik", "Zika", "Pvivax",
                "Pfalciparum"),
  name = "Co-exposure Profile",
  base_annotations = list(
    'Count' = intersection_size(
      mapping = aes(fill = after_stat(y))
    ) +
      scale_fill_gradient(
        low = "#8980BEFF",
        high = "#620A5DFF",
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
    upset_query(set = "Zika", fill = "#620A5DFF"),
    upset_query(set = "Chik", fill = "#620A5DFF"),
    upset_query(set = "Schistosomiasis", fill = "#8980BEFF"),
    upset_query(set = "Pfalciparum", fill = "#8980BEFF"),
    upset_query(set = "Pvivax", fill = "#8980BEFF")
  ),
  set_sizes = FALSE,
  matrix = intersection_matrix(
    geom = geom_point(size = 3),
    segment = geom_segment(size = 1)
  )
  # set_sizes = upset_set_size() +
  #   theme(
  #     panel.grid.major.x = element_blank(),
  #     panel.grid.minor.x = element_blank(),
  #     panel.grid.major.y = element_blank(),
  #     panel.grid.minor.y = element_blank()
  #   ) +
  #   ylab("Set")
) +
  ggtitle("Overall") +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

# Mostrar el gráfico
upset_plot 



# Upset Plot by age groups ------------------------------------------------

make_upset <- function(data, group_label) {
  ComplexUpset::upset(
    data,
    intersect = c(
      "Schistosomiasis",
      "Chik",
      "Zika",
      "Pvivax",
      "Pfalciparum"
    ),
    name = "Co-exposure Profiles",
    base_annotations = list(
      'Count' = intersection_size(mapping = aes(fill = after_stat(y)),
                                  text = list(size = 2.5)) +
        labs(y = "") +
        scale_fill_gradient(low = "#8980BEFF",
                            high = "#620A5DFF",
                            guide = "none") +
        theme(
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()
        )
    ),
    set_sizes = FALSE,
    matrix = intersection_matrix(
      geom = geom_point(size = 2.5),
      segment = geom_segment(size = 0.5)
    )
  ) + ggtitle(group_label) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      #axis.text.y = element_text(size = 6.5),
      axis.title.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    )
}

upset_by_age <- split(df_bin, df_bin$age_oms) %>% 
  imap(~ make_upset(.x, .y))

p1 <- upset_by_age[["Children and adolescent"]] +
  labs(x = "")

p2 <- upset_by_age[["Young adult"]] +
  labs (x = "")
p3 <- upset_by_age[["Middle-aged adult"]] 

p4 <- upset_by_age[["Older adult"]]


fig_2b <- wrap_plots(
  list(p1, p2, p3, p4),
  ncol = 2,
  nrow = 2,
  guides = "collect"
) 


fig2 <- wrap_plots(
  list(upset_plot, fig_2b),
  ncol = 2,
  nrow = 1
) 

fig2

# Guardar como imagen
ggsave("output/fig2_upset_plots.png",
       fig2,
       dpi = 600,
       width = 19,
       height = 10)


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
    Pvivax = if_else(
      pv_exposure == "Positive",
      1,
      0
    ),
    Pfalciparum = if_else(
      pf_exposure == "Positive",
      1,
      0
    )
  ) %>% 
  select(ffi_h_code, ffi_is_code,
         schistosomiasis, covid,
         chik, zika,
         Pvivax, Pfalciparum,
         ffi_gps_lat, ffi_gps_long) %>% 
  group_by(ffi_h_code) %>%
  summarise(
    prev_schisto  = sum(schistosomiasis, na.rm = TRUE) / n(),   
    prev_covid    = sum(covid, na.rm = TRUE) / n(),
    prev_chik     = sum(chik, na.rm = TRUE) / n(),
    prev_zika     = sum(zika, na.rm = TRUE) / n(),
    prev_pvivax  = sum(Pvivax, na.rm = TRUE) / n(),
    prev_pfalciparum = sum(Pfalciparum, na.rm = TRUE) / n(),
    lat  = first(ffi_gps_lat),
    long = first(ffi_gps_long)
  ) %>% 
  st_as_sf(
    coords = c("long", "lat"),
    crs = 4326
  )



# K = 1 ----------------------------------------------------------------------------------------
# Neighbours and spatial weights 
df_n1 <- df_spat %>% 
  mutate(
    nb = st_knn(geometry, k = 1),
    wt = st_weights(nb)
  )


# Spatial: Schistosomiasis
local_m_schisto1 <- df_n1 %>% 
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

local_m_schisto_rot1 <- local_m_schisto1 |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')
crop_coordinates_rot <- crop_coordinates |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')

bbox_crop1 <- st_bbox(c(xmin = -3179322 , xmax = 3199161,
                       ymin = -3696615 , ymax = -3603883),
                     crs = st_crs(local_m_schisto_rot1))

local_m_schisto_crop1 <- st_crop(local_m_schisto_rot1, bbox_crop1)

map_schisto1 <- local_m_schisto_crop1 %>%
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
  annotation_map_tile(type = "cartolight", zoom = 11) +
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

map_schisto1

# Spatial: P. falciparum
local_m_falciparum1 <- df_n1 %>% 
  mutate(
    lm_malaria = local_moran(prev_pfalciparum, nb, wt)
  ) %>% 
  unnest(lm_malaria)

local_m_malaria_rot1 <- local_m_falciparum1 |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')

bbox_crop1 <- st_bbox(c(xmin = -3179322 , xmax = 3199161,
                       ymin = -3696615 , ymax = -3603883),
                     crs = st_crs(local_m_malaria_rot1))

local_m_malaria_crop1 <- st_crop(local_m_malaria_rot1, bbox_crop1)

# Mapa

map_falciparum1 <- local_m_malaria_crop1 %>% 
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
  annotation_map_tile(type = "cartolight", zoom = 11) +
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

map_falciparum1

# Spatial: vivax
local_m_vivax1 <- df_n1 %>% 
  mutate(
    lm_vivax = local_moran(prev_pvivax, nb, wt)
  ) %>% 
  unnest(lm_vivax)

local_m_vivax_rot1 <- local_m_vivax1 |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')

bbox_crop1 <- st_bbox(c(xmin = -3179322 , xmax = 3199161,
                       ymin = -3696615 , ymax = -3603883),
                     crs = st_crs(local_m_vivax_rot1))

local_m_vivax_crop1 <- st_crop(local_m_vivax_rot1, bbox_crop1)

# Mapa

map_vivax1 <- local_m_vivax_crop1 %>% 
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
  annotation_map_tile(type = "cartolight", zoom = 11) +
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

map_vivax1

# Spatial: Zika
local_m_zika1 <- df_n1 %>% 
  mutate(
    lm_zika = local_moran(prev_zika, nb, wt)
  ) %>% 
  unnest(lm_zika)

local_m_zika_rot1 <- local_m_zika1 |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')

bbox_crop1 <- st_bbox(c(xmin = -3179322 , xmax = 3199161,
                       ymin = -3696615 , ymax = -3603883),
                     crs = st_crs(local_m_zika_rot1))

local_m_zika_crop1 <- st_crop(local_m_zika_rot1, bbox_crop1)

# Mapa

map_zika1 <- local_m_zika_crop1 %>% 
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
  annotation_map_tile(type = "cartolight", zoom = 11) +
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
  theme(
    legend.position = "none"
  ) +
  coord_sf(
    xlim = c(3172322, 3199761),
    ylim = c(-3698615, -3603883)
  )

map_zika1

# Spatial: Chik
local_m_chik1 <- df_n1 %>% 
  mutate(
    lm_chik = local_moran(prev_chik, nb, wt)
  ) %>% 
  unnest(lm_chik)

local_m_chik_rot1 <- local_m_chik1 |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')

bbox_crop1 <- st_bbox(c(xmin = -3179322 , xmax = 3199161,
                       ymin = -3696615 , ymax = -3603883),
                     crs = st_crs(local_m_chik_rot1))

local_m_chik_crop1 <- st_crop(local_m_chik_rot1, bbox_crop1)

# Mapa

map_chik1 <- local_m_chik_crop1 %>% 
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
  annotation_map_tile(type = "cartolight", zoom = 11) +
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

map_chik1

# Joining maps

map_five_join <- map_schisto1 + map_zika1 + map_falciparum1 + map_vivax1 + map_chik1 +
  plot_layout(
    ncol = 5
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

ggsave("output/five_maps_spatial_k1.png",
       map_five_join,
       dpi = 600,
       width = 15,
       height = 10)


# k = 2 ----------------------------------------------------------------------------------------

# Neighbours and spatial weights 
df_n2 <- df_spat %>% 
  mutate(
    nb = st_knn(geometry, k = 2),
    wt = st_weights(nb)
  )


# Spatial: Schistosomiasis
local_m_schisto2 <- df_n2 %>% 
  mutate(
    lm_schisto = local_moran(prev_schisto, nb, wt)
  ) %>% 
  unnest(lm_schisto)

# Mapa


local_m_schisto_rot2 <- local_m_schisto2 |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')

bbox_crop2 <- st_bbox(c(xmin = -3179322 , xmax = 3199161,
                       ymin = -3696615 , ymax = -3603883),
                     crs = st_crs(local_m_schisto_rot2))

local_m_schisto_crop2 <- st_crop(local_m_schisto_rot2, bbox_crop2)

map_schisto2 <- local_m_schisto_crop2 %>%
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
  annotation_map_tile(type = "cartolight", zoom = 11) +
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

map_schisto2

# Spatial: P. falciparum
local_m_falciparum2 <- df_n2 %>% 
  mutate(
    lm_malaria = local_moran(prev_pfalciparum, nb, wt)
  ) %>% 
  unnest(lm_malaria)

local_m_malaria_rot2 <- local_m_falciparum2 |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')

bbox_crop2 <- st_bbox(c(xmin = -3179322 , xmax = 3199161,
                       ymin = -3696615 , ymax = -3603883),
                     crs = st_crs(local_m_malaria_rot2))

local_m_malaria_crop2 <- st_crop(local_m_malaria_rot2, bbox_crop2)

# Mapa

map_falciparum2 <- local_m_malaria_crop2 %>% 
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
  annotation_map_tile(type = "cartolight", zoom = 11) +
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

map_falciparum2

# Spatial: vivax
local_m_vivax2 <- df_n2 %>% 
  mutate(
    lm_vivax = local_moran(prev_pvivax, nb, wt)
  ) %>% 
  unnest(lm_vivax)

local_m_vivax_rot2 <- local_m_vivax2 |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')

bbox_crop2 <- st_bbox(c(xmin = -3179322 , xmax = 3199161,
                       ymin = -3696615 , ymax = -3603883),
                     crs = st_crs(local_m_vivax_rot2))

local_m_vivax_crop2 <- st_crop(local_m_vivax_rot2, bbox_crop2)

# Mapa

map_vivax2 <- local_m_vivax_crop2 %>% 
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
  annotation_map_tile(type = "cartolight", zoom = 11) +
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

map_vivax2

# Spatial: Zika
local_m_zika2 <- df_n2 %>% 
  mutate(
    lm_zika = local_moran(prev_zika, nb, wt)
  ) %>% 
  unnest(lm_zika)

local_m_zika_rot2 <- local_m_zika2 |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')

bbox_crop2 <- st_bbox(c(xmin = -3179322 , xmax = 3199161,
                       ymin = -3696615 , ymax = -3603883),
                     crs = st_crs(local_m_zika_rot2))

local_m_zika_crop2 <- st_crop(local_m_zika_rot2, bbox_crop2)

# Mapa

map_zika2 <- local_m_zika_crop2 %>% 
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
  annotation_map_tile(type = "cartolight", zoom = 11) +
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
  theme(
    legend.position = "none"
  ) +
  coord_sf(
    xlim = c(3172322, 3199761),
    ylim = c(-3698615, -3603883)
  )

map_zika2

# Spatial: Chik
local_m_chik2 <- df_n2 %>% 
  mutate(
    lm_chik = local_moran(prev_chik, nb, wt)
  ) %>% 
  unnest(lm_chik)

local_m_chik_rot2 <- local_m_chik2 |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')

bbox_crop2 <- st_bbox(c(xmin = -3179322 , xmax = 3199161,
                       ymin = -3696615 , ymax = -3603883),
                     crs = st_crs(local_m_chik_rot2))

local_m_chik_crop2 <- st_crop(local_m_chik_rot2, bbox_crop2)

# Mapa

map_chik2 <- local_m_chik_crop2 %>% 
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
  annotation_map_tile(type = "cartolight", zoom = 11) +
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
    legend.position = "right"
  ) +
  coord_sf(
    xlim = c(3172322, 3199761),
    ylim = c(-3698615, -3603883)
  )

map_chik2

# Joining maps

map_five_join2 <- map_schisto2 + map_zika2 + map_falciparum2 + map_vivax2 + map_chik2 +
  plot_layout(
    ncol = 5
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

ggsave("output/five_maps_spatial_k2.png",
       map_five_join2,
       dpi = 600,
       width = 15,
       height = 10)



# k = 3 ----------------------------------------------------------------------------------------

# Neighbours and spatial weights 
df_n3 <- df_spat %>% 
  mutate(
    nb = st_knn(geometry, k = 3),
    wt = st_weights(nb)
  )


# Spatial: Schistosomiasis
local_m_schisto3 <- df_n3 %>% 
  mutate(
    lm_schisto = local_moran(prev_schisto, nb, wt)
  ) %>% 
  unnest(lm_schisto)

# Mapa


local_m_schisto_rot3 <- local_m_schisto3 |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')

bbox_crop3 <- st_bbox(c(xmin = -3179322 , xmax = 3199161,
                       ymin = -3696615 , ymax = -3603883),
                     crs = st_crs(local_m_schisto_rot3))

local_m_schisto_crop3 <- st_crop(local_m_schisto_rot3, bbox_crop3)

map_schisto3 <- local_m_schisto_crop3 %>%
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
  annotation_map_tile(type = "cartolight", zoom = 11) +
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

map_schisto3

# Spatial: P. falciparum
local_m_falciparum <- df_n %>% 
  mutate(
    lm_malaria = local_moran(prev_pfalciparum, nb, wt)
  ) %>% 
  unnest(lm_malaria)

local_m_malaria_rot <- local_m_falciparum |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')

bbox_crop <- st_bbox(c(xmin = -3179322 , xmax = 3199161,
                       ymin = -3696615 , ymax = -3603883),
                     crs = st_crs(local_m_malaria_rot))

local_m_malaria_crop <- st_crop(local_m_malaria_rot, bbox_crop)

# Mapa

map_falciparum <- local_m_malaria_crop %>% 
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
  annotation_map_tile(type = "cartolight", zoom = 11) +
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


# Spatial: vivax
local_m_vivax <- df_n %>% 
  mutate(
    lm_vivax = local_moran(prev_pvivax, nb, wt)
  ) %>% 
  unnest(lm_vivax)

local_m_vivax_rot <- local_m_vivax |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')

bbox_crop <- st_bbox(c(xmin = -3179322 , xmax = 3199161,
                       ymin = -3696615 , ymax = -3603883),
                     crs = st_crs(local_m_vivax_rot))

local_m_vivax_crop <- st_crop(local_m_vivax_rot, bbox_crop)

# Mapa

map_vivax <- local_m_vivax_crop %>% 
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
  annotation_map_tile(type = "cartolight", zoom = 11) +
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
  annotation_map_tile(type = "cartolight", zoom = 11) +
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
  annotation_map_tile(type = "cartolight", zoom = 11) +
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

map_five_join <- map_schisto + map_falciparum + map_vivax + map_chik + map_zika +
  plot_layout(
    ncol = 5
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

ggsave("output/five_maps_spatial.png",
       map_five_join,
       dpi = 600,
       width = 15,
       height = 10)



# K = 5 ----------------------------------------------------------------------------------------

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
  annotation_map_tile(type = "cartolight", zoom = 11) +
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


# Spatial: P. falciparum
local_m_falciparum <- df_n %>% 
  mutate(
    lm_malaria = local_moran(prev_pfalciparum, nb, wt)
  ) %>% 
  unnest(lm_malaria)

local_m_malaria_rot <- local_m_falciparum |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')

bbox_crop <- st_bbox(c(xmin = -3179322 , xmax = 3199161,
                       ymin = -3696615 , ymax = -3603883),
                     crs = st_crs(local_m_malaria_rot))

local_m_malaria_crop <- st_crop(local_m_malaria_rot, bbox_crop)

# Mapa

map_falciparum <- local_m_malaria_crop %>% 
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
  annotation_map_tile(type = "cartolight", zoom = 11) +
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


# Spatial: vivax
local_m_vivax <- df_n %>% 
  mutate(
    lm_vivax = local_moran(prev_pvivax, nb, wt)
  ) %>% 
  unnest(lm_vivax)

local_m_vivax_rot <- local_m_vivax |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')

bbox_crop <- st_bbox(c(xmin = -3179322 , xmax = 3199161,
                       ymin = -3696615 , ymax = -3603883),
                     crs = st_crs(local_m_vivax_rot))

local_m_vivax_crop <- st_crop(local_m_vivax_rot, bbox_crop)

# Mapa

map_vivax <- local_m_vivax_crop %>% 
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
  annotation_map_tile(type = "cartolight", zoom = 11) +
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
  annotation_map_tile(type = "cartolight", zoom = 11) +
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
  annotation_map_tile(type = "cartolight", zoom = 11) +
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

map_five_join <- map_schisto + map_falciparum + map_vivax + map_chik + map_zika +
  plot_layout(
    ncol = 5
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

ggsave("output/five_maps_spatial.png",
       map_five_join,
       dpi = 600,
       width = 15,
       height = 10)


