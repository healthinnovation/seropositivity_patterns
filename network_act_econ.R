library(tidyverse)
library(ComplexUpset)
library(sf)
library(sfdep)
library(ggspatial)
library(patchwork)
library(prettymapr)
library(igraph)
library(ggraph)

ffi <- read_csv("./data/mid/ffi_hh_ind_seropos.csv") %>% 
  mutate(
    act_econ = case_when(
      ffi_is_main_econ_act == "0" ~ "None",
      ffi_is_main_econ_act == "1" ~ "Day labourer",
      ffi_is_main_econ_act == "2" ~ "Wood extractor",
      ffi_is_main_econ_act == "3" ~ "Fisherman",
      ffi_is_main_econ_act == "4" ~ "Livestock farmer",
      ffi_is_main_econ_act == "5" ~ "Farmer",
      ffi_is_main_econ_act == "6" ~ "Trader",
      ffi_is_main_econ_act == "7" ~ "Housewife",
      ffi_is_main_econ_act == "8" ~ "Student",
      ffi_is_main_econ_act == "9" ~ "Motorcycle taxi driver",
      ffi_is_main_econ_act == "88" ~ "Other"
    )
  )
table(ffi$ffi_is_main_econ_act)
table(ffi$act_econ)
# Upset Plot --------------------------------------------------------------

#------------------------------------------------
# Data Preparation
df_bin <- ffi %>% 
  filter(act_econ=="Day labourer") %>% 
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
  pivot_longer(cols = all_of(enfermedades), names_to = "name", values_to = "count") %>% 
  filter(count>0)

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
    title = "A. Day labourer",
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

#------------------------------------------------
# Data Preparation
df_bin <- ffi %>% 
  filter(act_econ=="Wood extractor") %>% 
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
  pivot_longer(cols = all_of(enfermedades), names_to = "name", values_to = "count") %>% 
  filter(count>0)

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
    title = "B. Wood extractor",
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

#------------------------------------------------
# Data Preparation
df_bin <- ffi %>% 
  filter(act_econ=="Fisherman") %>% 
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
  pivot_longer(cols = all_of(enfermedades), names_to = "name", values_to = "count") %>% 
  filter(count>0)

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
c <- ggraph(grafo, layout = layout) +
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
    title = "C. Fisherman",
    size = "N° of seropositivity",
    edge_width = "N° of connections"
  ) +
  guides(
    size = guide_legend(
      override.aes = list(shape = 21, color = "white", fill = "#8980BEFF")
    )
  )

c

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

#------------------------------------------------
# Data Preparation
df_bin <- ffi %>% 
  filter(act_econ=="Livestock farmer") %>% 
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
  pivot_longer(cols = all_of(enfermedades), names_to = "name", values_to = "count") %>% 
  filter(count>0)

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
d <- ggraph(grafo, layout = layout) +
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
    title = "D. Livestock farmer",
    size = "N° of seropositivity",
    edge_width = "N° of connections"
  ) +
  guides(
    size = guide_legend(
      override.aes = list(shape = 21, color = "white", fill = "#8980BEFF")
    )
  )

d

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

#------------------------------------------------
# Data Preparation
df_bin <- ffi %>% 
  filter(act_econ=="Farmer") %>% 
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
  pivot_longer(cols = all_of(enfermedades), names_to = "name", values_to = "count") %>% 
  filter(count>0)

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
e <- ggraph(grafo, layout = layout) +
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
    title = "E. Farmer",
    size = "N° of seropositivity",
    edge_width = "N° of connections"
  ) +
  guides(
    size = guide_legend(
      override.aes = list(shape = 21, color = "white", fill = "#8980BEFF")
    )
  )

e

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

#------------------------------------------------
# Data Preparation
df_bin <- ffi %>% 
  filter(act_econ=="Trader") %>% 
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
  pivot_longer(cols = all_of(enfermedades), names_to = "name", values_to = "count") %>% 
  filter(count>0)

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
f <- ggraph(grafo, layout = layout) +
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
    title = "F. Trader",
    size = "N° of seropositivity",
    edge_width = "N° of connections"
  ) +
  guides(
    size = guide_legend(
      override.aes = list(shape = 21, color = "white", fill = "#8980BEFF")
    )
  )

f

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

#------------------------------------------------
# Data Preparation
df_bin <- ffi %>% 
  filter(act_econ=="Housewife") %>% 
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
  pivot_longer(cols = all_of(enfermedades), names_to = "name", values_to = "count") %>% 
  filter(count>0)

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
g <- ggraph(grafo, layout = layout) +
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
    title = "G. Housewife",
    size = "N° of seropositivity",
    edge_width = "N° of connections"
  ) +
  guides(
    size = guide_legend(
      override.aes = list(shape = 21, color = "white", fill = "#8980BEFF")
    )
  )

g

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

#------------------------------------------------
# Data Preparation
df_bin <- ffi %>% 
  filter(act_econ=="Student") %>% 
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
  pivot_longer(cols = all_of(enfermedades), names_to = "name", values_to = "count") %>% 
  filter(count>0)

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
h <- ggraph(grafo, layout = layout) +
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
    title = "H. Student",
    size = "N° of seropositivity",
    edge_width = "N° of connections"
  ) +
  guides(
    size = guide_legend(
      override.aes = list(shape = 21, color = "white", fill = "#8980BEFF")
    )
  )

h

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

#------------------------------------------------
# Data Preparation
df_bin <- ffi %>% 
  filter(act_econ=="Motorcycle taxi driver") %>% 
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
  pivot_longer(cols = all_of(enfermedades), names_to = "name", values_to = "count") %>% 
  filter(count>0)

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
i <- ggraph(grafo, layout = layout) +
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
    title = "I. Motorcycle taxi driver",
    size = "N° of seropositivity",
    edge_width = "N° of connections"
  ) +
  guides(
    size = guide_legend(
      override.aes = list(shape = 21, color = "white", fill = "#8980BEFF")
    )
  )

i

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

#------------------------------------------------
# Data Preparation
df_bin <- ffi %>% 
  filter(act_econ=="None") %>% 
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
  pivot_longer(cols = all_of(enfermedades), names_to = "name", values_to = "count") %>% 
  filter(count>0)

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
j <- ggraph(grafo, layout = layout) +
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
    title = "J. None",
    size = "N° of seropositivity",
    edge_width = "N° of connections"
  ) +
  guides(
    size = guide_legend(
      override.aes = list(shape = 21, color = "white", fill = "#8980BEFF")
    )
  )

j

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

#------------------------------------------------
# Data Preparation
df_bin <- ffi %>% 
  filter(act_econ=="Other") %>% 
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
  pivot_longer(cols = all_of(enfermedades), names_to = "name", values_to = "count") %>% 
  filter(count>0)

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
k <- ggraph(grafo, layout = layout) +
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
    title = "K. Other",
    size = "N° of seropositivity",
    edge_width = "N° of connections"
  ) +
  guides(
    size = guide_legend(
      override.aes = list(shape = 21, color = "white", fill = "#8980BEFF")
    )
  )

k

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
ggarrange(a,b,c,d,e,f,g,h,i,j,k, ncol = 4, nrow =3)
ggsave("output/network_act_econ.png", plot = last_plot(), 
       dpi = 600, height = 14.5, width = 21, bg = "white")
