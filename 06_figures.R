library(tidyverse)
library(ComplexUpset)
library(sf)
library(sfdep)
library(ggspatial)
library(patchwork)
library(prettymapr)
library(basemaps)
library(innovar)
library(scales)

ffi <- read_csv("./data/mid/ffi_hh_ind_seropos.csv") %>% 
  mutate(age_oms = case_when(
    age_cat %in% c("[0-10)", "[10-20)") ~ "Children and adolescent",
    age_cat %in% c("[20-30)", "[30-40)") ~ "Young adult",
    age_cat %in% c("[40-50)", "[50-60)") ~ "Middle-aged adult",
    age_cat %in% c("[60-70)", "[70+)")   ~ "Older adult",
    TRUE ~ NA_character_
    )
    )

# Spatial: Data Preparation 
ffi2 <- ffi %>% 
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
    ),
    coexposure = rowSums(
      across(
        c(
          schistosomiasis,
          covid,
          chik,
          zika,
          Pvivax,
          Pfalciparum)
        ), na.rm = T
      )
  ) %>% 
  select(ffi_h_code, ffi_is_code,
         schistosomiasis, covid,
         chik, zika,
         Pvivax, Pfalciparum, coexposure,
         ffi_gps_lat, ffi_gps_long) %>%
  group_by(ffi_h_code) %>%
  summarise(
    prev_schisto  = sum(schistosomiasis, na.rm = TRUE) / n(),   
    prev_covid    = sum(covid, na.rm = TRUE) / n(),
    prev_chik     = sum(chik, na.rm = TRUE) / n(),
    prev_zika     = sum(zika, na.rm = TRUE) / n(),
    prev_pvivax  = sum(Pvivax, na.rm = TRUE) / n(),
    prev_pfalciparum = sum(Pfalciparum, na.rm = TRUE) / n(),
    coexposure_tot = sum(coexposure, na.rm = TRUE),
    lat  = first(ffi_gps_lat),
    long = first(ffi_gps_long)
  ) 


df_spat <- ffi2 %>% 
  st_as_sf(
    coords = c("long", "lat"),
    crs = 4326
  )

# Figure 1 -------------------------------------------------------------------------------------


set_defaults(map_service = "carto",
             map_type = "light_no_labels")


# Changing projection
pts_3857 <- st_transform(df_spat, 3857) %>% 
  arrange(is.na(coexposure_tot), coexposure_tot)

# Adjusting buffer
bbox_3857 <- st_as_sfc(st_bbox(pts_3857), crs = 3857)
buffer_3857 <- st_buffer(bbox_3857, dist = 9000) 


# Adding adminsitrative limits
data("Peru", package = "innovar")

crop_coordinates <- Peru %>% 
  filter(dep == "LORETO",
         distr %in% c("INDIANA","BELEN")) %>% 
  st_transform(3857)

crop_coordinates_clip <- st_intersection(
  st_make_valid(crop_coordinates),
  st_make_valid(buffer_3857)
)

# >>> BBOX EXACTO DEL RECORTE para basemap y límites
bbox_clip <- st_bbox(crop_coordinates_clip)
bbox_geom <- st_as_sfc(bbox_clip, crs = 3857)

clip_union <- st_union(st_make_valid(crop_coordinates_clip))
mask_geom  <- st_make_valid(st_difference(bbox_geom, clip_union))


# (B) ***Borde verdadero***: línea del límite administrativo recortada al buffer
border_line_clip <- st_intersection(st_boundary(crop_coordinates), st_make_valid(buffer_3857))

# Map
exposures <- ggplot() +
  basemap_gglayer(bbox_geom) +
  coord_sf(expand = F) +
  scale_fill_identity() +
  geom_sf(data = mask_geom, fill = "white", colour = NA) +
  geom_sf(data = crop_coordinates_clip, #Para borde real -> border_line_clip
          fill = NA,
          color = "grey30",
          #linetype = 22,
          linewidth = 0.45) +
  geom_sf(data = pts_3857,
          aes(color = coexposure_tot),
          size = 4, alpha = 0.9) +
  scale_color_distiller(direction = 1,
                        palette = "Blues") +
  # annotation_north_arrow(location = "tr", style = north_arrow_nautical) +
  # coord_sf(expand = F) +
  # annotation_scale(location = "bl") +
  coord_sf(xlim = c(bbox_clip["xmin"], bbox_clip["xmax"]),
           ylim = c(bbox_clip["ymin"], bbox_clip["ymax"]),
           expand = FALSE) +
  theme_void(base_size = 11) +
  theme(
    legend.position = "none",
    panel.border = element_blank()
  )
  
exposures 



# Zona 1 --------------------------------------------------------------------------------

# === Centro y CÍRCULO en 3857 ===
center <- c(long = -72.90, lat = -3.50)

circle_3857 <- tibble(lon = center["long"], lat = center["lat"]) %>%
  st_as_sf(coords = c("lon","lat"), crs = 4326) %>%
  st_transform(3857) %>%
  st_buffer(dist = 15000)      # radio en METROS

# === Datos al MISMO CRS (3857) y recortes ===
pts_3857_ok <- st_transform(pts_3857, 3857)
borde_3857  <- st_transform(crop_coordinates_clip, 3857)

pts_in   <- st_filter(st_make_valid(pts_3857_ok), circle_3857)           # puntos dentro
borde_in <- st_intersection(st_make_valid(borde_3857), st_make_valid(circle_3857))

# === BBOX del zoom y MÁSCARA para recortar el basemap al círculo ===
bbox_zoom <- st_as_sfc(st_bbox(circle_3857), crs = 3857)
mask_zoom <- st_make_valid(st_difference(bbox_zoom, st_union(circle_3857)))
bb        <- st_bbox(circle_3857)

# Usa el mismo proveedor que en el mapa general
set_defaults(map_service = "carto", map_type = "light_no_labels")

# === Plot del zoom con basemap recortado ===
zoom1 <- ggplot() +
  basemap_gglayer(bbox_zoom) +
  scale_fill_identity() +
  geom_sf(data = mask_zoom, fill = "white", colour = NA) +
  geom_sf(data = borde_in, fill = NA, colour = "grey40",
          linetype = "22", linewidth = 0.45) +
  geom_sf(data = pts_in, aes(color = coexposure_tot),
          size = 3.5, alpha = 0.9) +
  geom_sf(data = st_boundary(circle_3857), colour = "black", linewidth = 0.4) +
  scale_color_distiller(direction = 1) +
  #scale_color_innova(palette = "npr", discrete = FALSE) +
  coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
           ylim = c(bb["ymin"], bb["ymax"]),
           expand = FALSE) +
  theme_void() +
  theme(legend.position = "none")

zoom1

# Zona 2 --------------------------------------------------------------------------------

# === Centro y CÍRCULO en 3857 ===
center2 <- c(long = -73.08, lat = -3.61)

circle_2 <- tibble(lon = center2["long"], lat = center2["lat"]) %>%
  st_as_sf(coords = c("lon","lat"), crs = 4326) %>%
  st_transform(3857) %>%
  st_buffer(dist = 4000)      # radio en METROS

# === Datos al MISMO CRS (3857) y recortes ===
pts_3857_ok2 <- st_transform(pts_3857, 3857)
borde_38572  <- st_transform(crop_coordinates_clip, 3857)

pts_in2   <- st_filter(st_make_valid(pts_3857_ok2), circle_2)           # puntos dentro
borde_in2 <- st_intersection(st_make_valid(borde_38572), st_make_valid(circle_2))

# === BBOX del zoom y MÁSCARA para recortar el basemap al círculo ===
bbox_zoom2 <- st_as_sfc(st_bbox(circle_2), crs = 3857)
mask_zoom2 <- st_make_valid(st_difference(bbox_zoom2, st_union(circle_2)))
bb2        <- st_bbox(circle_2)

# Usa el mismo proveedor que en el mapa general
set_defaults(map_service = "carto", map_type = "light_no_labels")

# === Plot del zoom con basemap recortado ===
zoom2 <- ggplot() +
  basemap_gglayer(bbox_zoom2) +
  scale_fill_identity() +
  geom_sf(data = mask_zoom2, fill = "white", colour = NA) +
  geom_sf(data = borde_in2, fill = NA, colour = "grey40",
          linetype = "22", linewidth = 0.45) +
  geom_sf(data = pts_in2, aes(color = coexposure_tot),
          size = 3.5, alpha = 0.9) +
  geom_sf(data = st_boundary(circle_2), colour = "black", linewidth = 0.4) +
  scale_color_distiller(direction = 1) +
  #scale_color_innova(palette = "npr", discrete = FALSE) +
  coord_sf(xlim = c(bb2["xmin"], bb2["xmax"]),
           ylim = c(bb2["ymin"], bb2["ymax"]),
           expand = FALSE) +
  theme_void() +
  theme(legend.position = "none")

zoom2

# Zona 3 --------------------------------------------------------------------------------

# === Centro y CÍRCULO en 3857 ===
center3 <- c(long = -73.25, lat = -3.77)

circle_3 <- tibble(lon = center3["long"], lat = center3["lat"]) %>%
  st_as_sf(coords = c("lon","lat"), crs = 4326) %>%
  st_transform(3857) %>%
  st_buffer(dist = 1500)      # radio en METROS

# === Datos al MISMO CRS (3857) y recortes ===
pts_3857_ok3 <- st_transform(pts_3857, 3857)
borde_38573  <- st_transform(crop_coordinates_clip, 3857)

pts_in3   <- st_filter(st_make_valid(pts_3857_ok3), circle_3)           # puntos dentro
borde_in3 <- st_intersection(st_make_valid(borde_38573), st_make_valid(circle_3))

# === BBOX del zoom y MÁSCARA para recortar el basemap al círculo ===
bbox_zoom3 <- st_as_sfc(st_bbox(circle_3), crs = 3857)
mask_zoom3 <- st_make_valid(st_difference(bbox_zoom3, st_union(circle_3)))
bb3        <- st_bbox(circle_3)

# Usa el mismo proveedor que en el mapa general
set_defaults(map_service = "carto", map_type = "light_no_labels")

# === Plot del zoom con basemap recortado ===
zoom3 <- ggplot() +
  basemap_gglayer(bbox_zoom3) +
  scale_fill_identity() +
  geom_sf(data = mask_zoom3, fill = "white", colour = NA) +
  geom_sf(data = borde_in3, fill = NA, colour = "grey40",
          linetype = "22", linewidth = 0.45) +
  geom_sf(data = pts_in3, aes(color = coexposure_tot),
          size = 3.5, alpha = 0.9) +
  geom_sf(data = st_boundary(circle_3), colour = "black", linewidth = 0.4) +
  scale_color_distiller(direction = 1) +
  #scale_color_innova(palette = "npr", discrete = FALSE) +
  coord_sf(xlim = c(bb3["xmin"], bb3["xmax"]),
           ylim = c(bb3["ymin"], bb3["ymax"]),
           expand = FALSE) +
  theme_void() +
  theme(legend.position = "none")

zoom3


# Zona 4 --------------------------------------------------------------------------------

# === Centro y CÍRCULO en 3857 ===
center4 <- c(long = -73.34, lat = -4.01)

circle_4 <- tibble(lon = center4["long"], lat = center4["lat"]) %>%
  st_as_sf(coords = c("lon","lat"), crs = 4326) %>%
  st_transform(3857) %>%
  st_buffer(dist = 10000)      # radio en METROS

# === Datos al MISMO CRS (3857) y recortes ===
pts_3857_ok4 <- st_transform(pts_3857, 3857)
borde_38574  <- st_transform(crop_coordinates_clip, 3857)

pts_in4   <- st_filter(st_make_valid(pts_3857_ok4), circle_4)           # puntos dentro
borde_in4 <- st_intersection(st_make_valid(borde_38574), st_make_valid(circle_4))

# === BBOX del zoom y MÁSCARA para recortar el basemap al círculo ===
bbox_zoom4 <- st_as_sfc(st_bbox(circle_4), crs = 3857)
mask_zoom4 <- st_make_valid(st_difference(bbox_zoom4, st_union(circle_4)))
bb4        <- st_bbox(circle_4)

# Usa el mismo proveedor que en el mapa general
set_defaults(map_service = "carto", map_type = "light_no_labels")

# === Plot del zoom con basemap recortado ===
zoom4 <- ggplot() +
  basemap_gglayer(bbox_zoom4) +
  scale_fill_identity() +
  geom_sf(data = mask_zoom4, fill = "white", colour = NA) +
  geom_sf(data = borde_in4, fill = NA, colour = "grey40",
          linetype = "22", linewidth = 0.45) +
  geom_sf(data = pts_in4, aes(color = coexposure_tot),
          size = 3.5, alpha = 0.9) +
  geom_sf(data = st_boundary(circle_4), colour = "black", linewidth = 0.4) +
  scale_color_distiller(direction = 1) +
  #scale_color_innova(palette = "npr", discrete = FALSE) +
  coord_sf(xlim = c(bb4["xmin"], bb4["xmax"]),
           ylim = c(bb4["ymin"], bb4["ymax"]),
           expand = FALSE) +
  theme_void() +
  theme(legend.position = "none")

zoom4



# Joining Maps ---------------------------------------------------------------------------------

base <- exposures + theme(plot.margin = margin(t = 60, r = 260, b = 60, l = 240))

fig1 <- base +
  inset_element(zoom1,
                left = 1.3, right = 1.8,
                bottom = 0.75, top = 1.09,
                align_to = "plot", clip = FALSE
  ) +
  inset_element(zoom2,
                left = 1.05, right = 1.55,
                bottom = 0.4, top = 0.74,
                align_to = "plot", clip = FALSE
  ) +
  inset_element(zoom3,
                left = -0.32, right = 0.12,
                bottom = 0.6, top = 0.94,
                align_to = "plot", clip = FALSE
  ) +
  inset_element(zoom4,
                left = -0.62, right = -0.12,
                bottom = 0.08, top = 0.42,
                align_to = "plot", clip = FALSE
  )


exposures_leg <- exposures +
  scale_color_distiller(
    direction = 1,
    palette = "Blues",
    name = "Number of \nco-exposures"   # título de la leyenda
  ) +
  guides(
    colour = guide_colourbar(     # dale formato a la barra
      ticks = FALSE,
      barheight = unit(80, "pt"),
      barwidth  = unit(6,  "pt")
    )
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 9)
  )


# 2) Extraer la leyenda como grob y convertirla en un plot
leg_grob <- cowplot::get_legend(exposures_leg)
leg_plot <- cowplot::ggdraw(leg_grob)


fig1_with_leg <- fig1 +
  inset_element(
    leg_plot,
    left = 1.05, right = 1.55,   
    bottom = 0.06, top = 0.38,   
    align_to = "plot",
    clip = FALSE
  )



library(grid)
final <- cowplot::ggdraw(fig1_with_leg) +
  geom_curve(
    aes(x = 0.61, y = 0.73, xend = 0.76, yend = 0.8),
    curvature = 0.2,
    arrow = arrow(type = "closed", length = unit(3, "mm")),
    size = 0.6, colour = "black"
  ) +
  geom_curve(
    aes(x = 0.50, y = 0.66, xend = 0.68, yend = 0.55),
    curvature = 0.2,
    arrow = arrow(type = "closed", length = unit(3, "mm")),
    size = 0.6, colour = "black"
  ) +
  geom_curve(
    aes(x = 0.40, y = 0.53, xend = 0.36, yend = 0.61),
    curvature = 0.2,
    arrow = arrow(type = "closed", length = unit(3, "mm")),
    size = 0.6, colour = "black"
  ) +
  geom_curve(
    aes(x = 0.33, y = 0.31, xend = 0.28, yend = 0.35),
    curvature = 0.2,
    arrow = arrow(type = "closed", length = unit(3, "mm")),
    size = 0.6, colour = "black"
  )
final 

ggsave(
  filename = "./output/figure1.png",
  plot = final,
  width = 13,
  height = 7,
  dpi = 1200
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
  min_size = 1,
  min_degree = 1,
  base_annotations = list(
    'Count' = intersection_size(
      mapping = aes(fill = after_stat(y))
    ) +
      scale_fill_gradient(
        low = "#78B1D3",
        high = "#3F8BBA",
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
    min_degree = 1,
    base_annotations = list(
      'Count' = intersection_size(mapping = aes(fill = after_stat(y)),
                                  text = list(size = 2.5)) +
        labs(y = "") +
        scale_fill_gradient(low = "#78B1D3",
                            high = "#3F8BBA",
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
    legend.position = "right"
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
  theme(
    legend.position = "none"
  ) +
  coord_sf(
    xlim = c(3172322, 3199761),
    ylim = c(-3698615, -3603883)
  )

map_chik1

# Joining maps

map_five_join <- map_schisto1 + map_falciparum1 + map_vivax1 + map_chik1 + map_zika1 +
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
    legend.position = "right"
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
    legend.position = "none"
  ) +
  coord_sf(
    xlim = c(3172322, 3199761),
    ylim = c(-3698615, -3603883)
  )

map_chik2

# Joining maps

map_five_join2 <- map_schisto2 + map_falciparum2 + map_vivax2 + map_chik2 + map_zika2 +
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
local_m_falciparum3 <- df_n3 %>% 
  mutate(
    lm_malaria = local_moran(prev_pfalciparum, nb, wt)
  ) %>% 
  unnest(lm_malaria)

local_m_malaria_rot3 <- local_m_falciparum3 |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')

bbox_crop3 <- st_bbox(c(xmin = -3179322 , xmax = 3199161,
                       ymin = -3696615 , ymax = -3603883),
                     crs = st_crs(local_m_malaria_rot3))

local_m_malaria_crop3 <- st_crop(local_m_malaria_rot3, bbox_crop3)

# Mapa

map_falciparum3 <- local_m_malaria_crop3 %>% 
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

map_falciparum3

# Spatial: vivax
local_m_vivax3 <- df_n3 %>% 
  mutate(
    lm_vivax = local_moran(prev_pvivax, nb, wt)
  ) %>% 
  unnest(lm_vivax)

local_m_vivax_rot3 <- local_m_vivax3 |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')

bbox_crop3 <- st_bbox(c(xmin = -3179322 , xmax = 3199161,
                       ymin = -3696615 , ymax = -3603883),
                     crs = st_crs(local_m_vivax_rot3))

local_m_vivax_crop3 <- st_crop(local_m_vivax_rot3, bbox_crop3)

# Mapa

map_vivax3 <- local_m_vivax_crop3 %>% 
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

map_vivax3

# Spatial: Zika
local_m_zika3 <- df_n3 %>% 
  mutate(
    lm_zika = local_moran(prev_zika, nb, wt)
  ) %>% 
  unnest(lm_zika)

local_m_zika_rot3 <- local_m_zika3 |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')

bbox_crop3 <- st_bbox(c(xmin = -3179322 , xmax = 3199161,
                       ymin = -3696615 , ymax = -3603883),
                     crs = st_crs(local_m_zika_rot3))

local_m_zika_crop3 <- st_crop(local_m_zika_rot3, bbox_crop3)

# Mapa

map_zika3 <- local_m_zika_crop3 %>% 
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

map_zika3

# Spatial: Chik
local_m_chik3 <- df_n3 %>% 
  mutate(
    lm_chik = local_moran(prev_chik, nb, wt)
  ) %>% 
  unnest(lm_chik)

local_m_chik_rot3 <- local_m_chik3 |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')

bbox_crop3 <- st_bbox(c(xmin = -3179322 , xmax = 3199161,
                       ymin = -3696615 , ymax = -3603883),
                     crs = st_crs(local_m_chik_rot3))

local_m_chik_crop3 <- st_crop(local_m_chik_rot3, bbox_crop3)

# Mapa

map_chik3 <- local_m_chik_crop3 %>% 
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

map_chik3

# Joining maps

map_five_join3 <- map_schisto3 + map_falciparum3 + map_vivax3 + map_chik3 + map_zika3 +
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

ggsave("output/five_maps_spatial_k3.png",
       map_five_join3,
       dpi = 600,
       width = 15,
       height = 10)



# k = 4 -------------------------------------------------------------------

# Neighbours and spatial weights 
df_n4 <- df_spat %>% 
  mutate(
    nb = st_knn(geometry, k = 4),
    wt = st_weights(nb)
  )


# Spatial: Schistosomiasis
local_m_schisto4 <- df_n4 %>% 
  mutate(
    lm_schisto = local_moran(prev_schisto, nb, wt)
  ) %>% 
  unnest(lm_schisto)

# Mapa

local_m_schisto_rot4 <- local_m_schisto4 |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')

bbox_crop4 <- st_bbox(c(xmin = -3179322 , xmax = 3199161,
                       ymin = -3696615 , ymax = -3603883),
                     crs = st_crs(local_m_schisto_rot4))

local_m_schisto_crop4 <- st_crop(local_m_schisto_rot4, bbox_crop4)

map_schisto4 <- local_m_schisto_crop4 %>%
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

map_schisto4

# Spatial: P. falciparum
local_m_falciparum4 <- df_n4 %>% 
  mutate(
    lm_malaria = local_moran(prev_pfalciparum, nb, wt)
  ) %>% 
  unnest(lm_malaria)

local_m_malaria_rot4 <- local_m_falciparum4 |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')

bbox_crop4 <- st_bbox(c(xmin = -3179322 , xmax = 3199161,
                       ymin = -3696615 , ymax = -3603883),
                     crs = st_crs(local_m_malaria_rot4))

local_m_malaria_crop4 <- st_crop(local_m_malaria_rot4, bbox_crop4)

# Mapa

map_falciparum4 <- local_m_malaria_crop4 %>% 
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

map_falciparum4

# Spatial: vivax
local_m_vivax4 <- df_n4 %>% 
  mutate(
    lm_vivax = local_moran(prev_pvivax, nb, wt)
  ) %>% 
  unnest(lm_vivax)

local_m_vivax_rot4 <- local_m_vivax4 |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')

bbox_crop4 <- st_bbox(c(xmin = -3179322 , xmax = 3199161,
                       ymin = -3696615 , ymax = -3603883),
                     crs = st_crs(local_m_vivax_rot4))

local_m_vivax_crop4 <- st_crop(local_m_vivax_rot4, bbox_crop4)

# Mapa

map_vivax4 <- local_m_vivax_crop4 %>% 
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

map_vivax4

# Spatial: Zika
local_m_zika4 <- df_n4 %>% 
  mutate(
    lm_zika = local_moran(prev_zika, nb, wt)
  ) %>% 
  unnest(lm_zika)

local_m_zika_rot4 <- local_m_zika4 |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')

bbox_crop4 <- st_bbox(c(xmin = -3179322 , xmax = 3199161,
                       ymin = -3696615 , ymax = -3603883),
                     crs = st_crs(local_m_zika_rot4))

local_m_zika_crop4 <- st_crop(local_m_zika_rot4, bbox_crop4)

# Mapa

map_zika4 <- local_m_zika_crop4 %>% 
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

map_zika4

# Spatial: Chik
local_m_chik4 <- df_n4 %>% 
  mutate(
    lm_chik = local_moran(prev_chik, nb, wt)
  ) %>% 
  unnest(lm_chik)

local_m_chik_rot4 <- local_m_chik4 |> 
  st_transform('+proj=omerc +lat_0=40 +lonc=-74 +gamma=-40')

bbox_crop4 <- st_bbox(c(xmin = -3179322 , xmax = 3199161,
                       ymin = -3696615 , ymax = -3603883),
                     crs = st_crs(local_m_chik_rot4))

local_m_chik_crop4 <- st_crop(local_m_chik_rot4, bbox_crop4)

# Mapa

map_chik4 <- local_m_chik_crop4 %>% 
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

map_chik4

# Joining maps

map_five_join4 <- map_schisto4 + map_falciparum4 + map_vivax4 + map_chik4 + map_zika4 +
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

ggsave("output/five_maps_spatial_k4.png",
       map_five_join4,
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


