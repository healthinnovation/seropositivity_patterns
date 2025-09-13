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

# Upset Plot --------------------------------------------------------------

# Data Preparation
df_bin <- ffi %>% 
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
  separate(pair, into = c("from", "to"), sep = "_")

grafo <- graph_from_data_frame(d = enlaces, vertices = nodos, directed = FALSE)

ggraph(grafo, layout = "circle") +
  geom_edge_link(aes(width = weight), color = "grey50") +
  geom_node_point(aes(size = count), color = "skyblue") +
  geom_node_text(aes(label = sub("_(Female|Male)$", "", name)), repel = TRUE) +
  scale_edge_width(range = c(0.5, 3)) +
  scale_size(range = c(5, 15)) +
  theme_void()

library(here)
library(tidyverse)
library(igraph)
library(igraphdata)

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

centrality_summary

# Calculate edge betweenness
edge_betweenness(g)

# Extract edges with the highest betweenness centrality
eb <- edge_betweenness(g)
as_edgelist(g)[which.max(eb), ]

# Obtaining the line graph G' of a graph G
lg <- make_line_graph(g)
lg

#Network cohesion
# Count frequency of cliques of each size
table(lengths(cliques(g)))

# Count the frequency of maximal cliques of each size
table(lengths(max_cliques(g)))

# Calculate density of the entire graph
edge_density(g)

# Calculate clustering coefficient of the graph
transitivity(g)

# Calculate reciprocity of the `enron` network
reciprocity(g)

# Check if the network is connected
is_connected(g)

# Decompose the graph and count number of nodes in each component
comps <- decompose(g)
table(lengths(comps))

# Extract giant component
yeast_gc <- decompose(g)[[1]]

# Get average path length and diameter
mean_distance(yeast_gc)

diameter(yeast_gc)

# Get clustering coefficient
transitivity(yeast_gc)

# Get articulation points
ap <- articulation_points(yeast_gc)
length(ap)

# Detect communities with Infomap
cl <- cluster_infomap(g)
cl

# Vector of cluster membership for each node
membership(cl)

# List of nodes that belong to each cluster
communities(cl)

# Visualize graph with communities highlighted
plot(cl, g)

# Detect community with Louvain
cl2 <- cluster_louvain(g)

# Compare community structures
compare(cl, cl2)

plot(cl2, g)

