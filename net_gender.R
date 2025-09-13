library(tidyverse)
library(igraph)
library(ggraph)
library(ggpubr)

ffi <- read_csv("./data/mid/ffi_hh_ind_seropos.csv")

# --- Función auxiliar para preparar grafo por subconjunto ----
prep_graph <- function(df){
  enfermedades <- c("Schistosomiasis","Zika","Chik","Malariavivax","Malariafalciparum")
  
  nodos <- df %>%
    summarise(across(all_of(enfermedades), ~sum(. == "Si")), .groups = "drop") %>%
    pivot_longer(cols = all_of(enfermedades), names_to = "name", values_to = "count")
  
  enlaces <- df %>%
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
      cols = starts_with(c("Schistosomiasis","Zika","Chik","Malariavivax","Malariafalciparum")),
      names_to = "pair", values_to = "weight"
    ) %>%
    separate(pair, into = c("from","to"), sep = "_") %>%
    filter(weight > 0)
  
  g <- graph_from_data_frame(d = enlaces, vertices = nodos, directed = FALSE)
  
  # atributos de vértices
  labels_custom <- c(
    "Schistosomiasis"   = "Schistosomiasis",
    "Zika"              = "Zika",
    "Chik"              = "Chikungunya",
    "Malariavivax"      = "italic('P. vivax')",
    "Malariafalciparum" = "italic('P. falciparum')"
  )
  V(g)$label_custom <- labels_custom[V(g)$name]
  
  # resaltar el mayor eigenvector
  ev <- eigen_centrality(g)$vector
  top_node <- names(ev)[which.max(ev)]
  V(g)$grupo <- ifelse(V(g)$name == top_node, "Menor closeness", "Otros nodos")
  
  # layout con atracción por weight (más peso = más cercanía)
  lay <- layout_with_fr(g, weights = E(g)$weight)
  
  list(g = g, layout = lay)
}

# --- Datos binarios para cada género ----
mk_bin <- function(sex){
  ffi %>%
    filter(gender == sex) %>%
    mutate(
      Schistosomiasis = if_else(SEA_pos == 1 | X229.E.NP_pos == 1, "Si","No"),
      Chik            = if_else(Chik.E1_pos == 1, "Si","No"),
      Zika            = if_else(Zika.NS1_pos == 1, "Si","No"),
      Malariavivax    = if_else(pv_exposure == "Positive", "Si","No"),
      Malariafalciparum = if_else(pf_exposure == "Positive", "Si","No")
    ) %>%
    select(ffi_is_code, Schistosomiasis, Chik, Zika, Malariavivax, Malariafalciparum) %>%
    drop_na()
}

df_f <- mk_bin("Female")
df_m <- mk_bin("Male")

pg_f <- prep_graph(df_f)
pg_m <- prep_graph(df_m)

g_f <- pg_f$g; g_m <- pg_m$g
lay_f <- pg_f$layout; lay_m <- pg_m$layout

# --- RANGOS GLOBALES (claves para una única leyenda proporcional) ----
global_edge_limits  <- range(c(E(g_f)$weight,  E(g_m)$weight),  na.rm = TRUE)
global_node_limits  <- range(c(V(g_f)$count,   V(g_m)$count),   na.rm = TRUE)

# --- ESCALAS GLOBALES (compartidas por a y b) ----
edge_width_global <- scale_edge_width(
  range  = c(0.5, 3),
  limits = global_edge_limits,
  name   = "N° of connections"
)

node_size_global <- scale_size_continuous(
  range  = c(5, 15),
  limits = global_node_limits,
  name   = "N° of seropositivity"
)

fill_nodes_global <- scale_fill_manual(
  values = c("Menor closeness" = "#620A5DFF", "Otros nodos" = "#8980BEFF"),
  guide  = "none"  # sin leyenda de color para nodos
)

# --- Función para dibujar usando las escalas globales ----
plot_graph <- function(g, lay, title_text){
  ggraph(g, layout = lay) +
    # aristas: grosor por weight, color CONSTANTE
    geom_edge_link(aes(width = weight), color = "#A9C1DDFF") +
    # nodos: tamaño por count, paleta de 2 colores (máx eigenvector vs resto)
    geom_node_point(aes(size = count, fill = grupo), shape = 21, color = "white") +
    geom_node_text(
      aes(label = label_custom),
      repel = TRUE, parse = TRUE,
      box.padding = unit(0.6, "lines"),
      point.padding = unit(0.5, "lines"),
      segment.color = NA
    ) +
    edge_width_global +
    node_size_global +
    fill_nodes_global +
    theme_void() +
    labs(title = title_text) +
    # asegurar que cada panel exponga las mismas guías
    guides(
      width = guide_legend(order = 1),  # edge_width
      size  = guide_legend(
        order = 2,
        override.aes = list(shape = 21, color = "white", fill = "#8980BEFF")
      )
    ) +
    theme(legend.position = "right")
}

a <- plot_graph(g_f, lay_f, "A. Female")
b <- plot_graph(g_m, lay_m, "B. Male")

# --- Unir con UNA sola leyenda global (porque comparten escalas y guías) ----
p <- ggarrange(a, b, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
p

ggsave("output/network_gender.png", plot = p, dpi = 600, height = 5, width = 11, bg = "white")
