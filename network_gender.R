library(tidyverse)
library(ComplexUpset)
library(sf)
library(sfdep)
library(ggspatial)
library(patchwork)
library(prettymapr)
library(igraph)
library(ggraph)

ffi <- read_csv("./data/mid/ffi_hh_ind_seropos.csv") 
table(ffi$ffi_is_main_econ_act)
# Upset Plot --------------------------------------------------------------

# Data Preparation
df_bin <- ffi %>% 
  filter(gender=="Female") %>% 
  mutate(
    Schistosomiasis = if_else(
      SEA_pos == 1 | X229.E.NP_pos == 1,
      "Si",
      "No"
    ),
    Chik = if_else(
      Chik.E1_pos == 1, 
      "Si", 
      "No"
    ),
    Zika = if_else(
      Zika.NS1_pos == 1,
      "Si",
      "No"
    ),
    Malariavivax = if_else(
      pv_exposure == "Positive",
      "Si",
      "No"
    ),
    Malariafalciparum = if_else(
      pf_exposure == "Positive",
      "Si",
      "No"
    )
  )%>% 
  select(ffi_is_code, Schistosomiasis, Chik, Zika,
         Malariavivax, Malariafalciparum) %>% 
  drop_na()

enfermedades <- c("Schistosomiasis", "Zika", "Chik", "Malariavivax", "Malariafalciparum")

# Calcular nodos y aristas por sexo
nodos <- df_bin %>%
  summarise(across(all_of(enfermedades), ~sum(. == "Si")), .groups = "drop") %>%
  pivot_longer(cols = all_of(enfermedades), names_to = "name", values_to = "count")

enlaces <- df_bin %>%
  summarise(
    Schistosomiasis_Zika = sum(Schistosomiasis == "Si" & Zika == "Si"),
    Schistosomiasis_Chik = sum(Schistosomiasis == "Si" & Chik == "Si"),
    Schistosomiasis_Malariavivax = sum(Schistosomiasis == "Si" & Malariavivax == "Si"),
    Schistosomiasis_Malariafalciparum = sum(Schistosomiasis == "Si" & Malariafalciparum == "Si"),
    Zika_Chik = sum(Zika == "Si" & Chik == "Si"),
    Zika_Malariavivax = sum(Zika == "Si" & Malariavivax == "Si"),
    Zika_Malariafalciparum = sum(Zika == "Si" & Malariafalciparum == "Si"),
    Chik_Malariavivax = sum(Chik == "Si" & Malariavivax == "Si"),
    Chik_Malariafalciparum = sum(Chik == "Si" & Malariafalciparum == "Si"),
    Malariavivax_Malariafalciparum = sum(Malariavivax == "Si" & Malariafalciparum == "Si"),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = starts_with(c("Schistosomiasis", "Zika", "Chik", "Malariavivax", "Malariafalciparum")),
    names_to = "pair", values_to = "weight"
  ) %>%
  separate(pair, into = c("from", "to"), sep = "_") %>% 
  filter(weight>0)

grafo <- graph_from_data_frame(d = enlaces, vertices = nodos, directed = FALSE)

E(grafo)$weight_dist <- 1 / E(grafo)$weight

# Usar layout con distancia inversa al peso (más peso = más cercanía)
layout <- layout_with_fr(grafo, weights = E(grafo)$weight)  # se usa "weights" para atracción
# O también layout_with_kk(grafo, weights = E(grafo)$weight)

# Graficar con ggraph
closeness_vals <- eigen_centrality(grafo)$vector

# Identificar nodo con menor closeness
min_node <- names(closeness_vals)[which.max(closeness_vals)]

# Labels customizados en formato de texto con italics solo en los Plasmodium
labels_custom <- c(
  "Schistosomiasis"   = "Schistosomiasis",
  "Zika"              = "Zika",
  "Chik"              = "Chikungunya",
  "Malariavivax"      = "italic('P. vivax')",
  "Malariafalciparum" = "italic('P. falciparum')"
)

# Guardamos como texto (character)
V(grafo)$label_custom <- labels_custom[V(grafo)$name]
V(grafo)$grupo <- ifelse(V(grafo)$name == min_node, 
                         "Menor closeness", 
                         "Otros nodos")

# Graficar con parse=TRUE para interpretar los italics
a <- ggraph(grafo, layout = layout) +
  geom_edge_link(aes(width = weight), color = "#A9C1DDFF") +
  geom_node_point(aes(size = count, fill = grupo), shape = 21, color = "white") +
  geom_node_text(
    aes(label = label_custom),
    repel = TRUE,
    parse = TRUE,
    box.padding = unit(0.6, "lines"),   
    point.padding = unit(0.5, "lines"),
    segment.color = NA   # <- elimina las líneas de unión
  ) +
  scale_edge_width(range = c(0.5, 3)) +
  scale_size(range = c(5, 15), guide = "legend") +
  scale_fill_manual(values = c("Menor closeness" = "#620A5DFF", 
                               "Otros nodos" = "#8980BEFF"), guide = "none") +
  theme_void() +
  labs(
    title = "A. Female",
    size = "N° of seropositivity",
    edge_width = "N° of connections"
  ) +
  guides(
    size = guide_legend(
      override.aes = list(shape = 21, color = "white", fill = "#8980BEFF")
    )
  )

a

library(htmlTable)

g <- grafo

# Degree
c1 <- degree(g)

# Closeness
c2 <- closeness(g)

# Betweenness
c3 <- betweenness(g)

# Eigenvector centrality
c4 <- eigen_centrality(g)$vector

# Harmonic centrality
c5 <- harmonic_centrality(g)

# Show all measures
centrality_summary <- data.frame(
  node = as.character(V(g)),
  degree = c1,
  closeness = c2,
  betweenness = c3,
  eigenvector = c4,
  harmonic = c5
)

htmlTable(centrality_summary)

#----------------------------------------
df_bin <- ffi %>% 
  filter(gender=="Male") %>% 
  mutate(
    Schistosomiasis = if_else(
      SEA_pos == 1 | X229.E.NP_pos == 1,
      "Si",
      "No"
    ),
    Chik = if_else(
      Chik.E1_pos == 1, 
      "Si", 
      "No"
    ),
    Zika = if_else(
      Zika.NS1_pos == 1,
      "Si",
      "No"
    ),
    Malariavivax = if_else(
      pv_exposure == "Positive",
      "Si",
      "No"
    ),
    Malariafalciparum = if_else(
      pf_exposure == "Positive",
      "Si",
      "No"
    )
  )%>% 
  select(ffi_is_code, Schistosomiasis, Chik, Zika,
         Malariavivax, Malariafalciparum) %>% 
  drop_na()

enfermedades <- c("Schistosomiasis", "Zika", "Chik", "Malariavivax", "Malariafalciparum")

# Calcular nodos y aristas por sexo
nodos <- df_bin %>%
  summarise(across(all_of(enfermedades), ~sum(. == "Si")), .groups = "drop") %>%
  pivot_longer(cols = all_of(enfermedades), names_to = "name", values_to = "count")

enlaces <- df_bin %>%
  summarise(
    Schistosomiasis_Zika = sum(Schistosomiasis == "Si" & Zika == "Si"),
    Schistosomiasis_Chik = sum(Schistosomiasis == "Si" & Chik == "Si"),
    Schistosomiasis_Malariavivax = sum(Schistosomiasis == "Si" & Malariavivax == "Si"),
    Schistosomiasis_Malariafalciparum = sum(Schistosomiasis == "Si" & Malariafalciparum == "Si"),
    Zika_Chik = sum(Zika == "Si" & Chik == "Si"),
    Zika_Malariavivax = sum(Zika == "Si" & Malariavivax == "Si"),
    Zika_Malariafalciparum = sum(Zika == "Si" & Malariafalciparum == "Si"),
    Chik_Malariavivax = sum(Chik == "Si" & Malariavivax == "Si"),
    Chik_Malariafalciparum = sum(Chik == "Si" & Malariafalciparum == "Si"),
    Malariavivax_Malariafalciparum = sum(Malariavivax == "Si" & Malariafalciparum == "Si"),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = starts_with(c("Schistosomiasis", "Zika", "Chik", "Malariavivax", "Malariafalciparum")),
    names_to = "pair", values_to = "weight"
  ) %>%
  separate(pair, into = c("from", "to"), sep = "_") %>% 
  filter(weight>0)

grafo <- graph_from_data_frame(d = enlaces, vertices = nodos, directed = FALSE)

E(grafo)$weight_dist <- 1 / E(grafo)$weight

# Usar layout con distancia inversa al peso (más peso = más cercanía)
layout <- layout_with_fr(grafo, weights = E(grafo)$weight)  # se usa "weights" para atracción
# O también layout_with_kk(grafo, weights = E(grafo)$weight)

# Graficar con ggraph
closeness_vals <- eigen_centrality(grafo)$vector

# Identificar nodo con menor closeness
min_node <- names(closeness_vals)[which.max(closeness_vals)]

# Labels customizados en formato de texto con italics solo en los Plasmodium
labels_custom <- c(
  "Schistosomiasis"   = "Schistosomiasis",
  "Zika"              = "Zika",
  "Chik"              = "Chikungunya",
  "Malariavivax"      = "italic('P. vivax')",
  "Malariafalciparum" = "italic('P. falciparum')"
)

# Guardamos como texto (character)
V(grafo)$label_custom <- labels_custom[V(grafo)$name]
V(grafo)$grupo <- ifelse(V(grafo)$name == min_node, 
                         "Menor closeness", 
                         "Otros nodos")

# Graficar con parse=TRUE para interpretar los italics
b <- ggraph(grafo, layout = layout) +
  geom_edge_link(aes(width = weight), color = "#A9C1DDFF") +
  geom_node_point(aes(size = count, fill = grupo), shape = 21, color = "white") +
  geom_node_text(
    aes(label = label_custom),
    repel = TRUE,
    parse = TRUE,
    box.padding = unit(0.6, "lines"),   
    point.padding = unit(0.5, "lines"),
    segment.color = NA   # <- elimina las líneas de unión
  ) +
  scale_edge_width(range = c(0.5, 3)) +
  scale_size(range = c(5, 15), guide = "legend") +
  scale_fill_manual(values = c("Menor closeness" = "#620A5DFF", 
                               "Otros nodos" = "#8980BEFF"), guide = "none") +
  theme_void() +
  labs(
    title = "B. Male",
    size = "N° of seropositivity",
    edge_width = "N° of connections"
  ) +
  guides(
    size = guide_legend(
      override.aes = list(shape = 21, color = "white", fill = "#8980BEFF")
    )
  )

b

library(htmlTable)

g <- grafo

# Degree
c1 <- degree(g)

# Closeness
c2 <- closeness(g)

# Betweenness
c3 <- betweenness(g)

# Eigenvector centrality
c4 <- eigen_centrality(g)$vector

# Harmonic centrality
c5 <- harmonic_centrality(g)

# Show all measures
centrality_summary <- data.frame(
  node = as.character(V(g)),
  degree = c1,
  closeness = c2,
  betweenness = c3,
  eigenvector = c4,
  harmonic = c5
)

htmlTable(centrality_summary)

library(ggpubr)
ggarrange(a,b, ncol = 2, nrow =1)
ggsave("output/network_gender.png", plot = last_plot(), 
       dpi = 600, height = 5, width = 12, bg = "white")
