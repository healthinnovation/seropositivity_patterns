pacman::p_load("sf","spatstat","dplyr")

households <- st_read("points.gpkg") |> 
  st_transform(crs = 32718)

# Coordinate extracting
coords <- st_coordinates(households)
x <- coords[, 1]
y <- coords[, 2]

# New bounding box for spatstat
bbox <- st_bbox(households)
ventana <- owin(xrange = c(bbox["xmin"], bbox["xmax"]),
                yrange = c(bbox["ymin"], bbox["ymax"]))

# Create a point pattern object
puntos <- ppp(x, y, window = ventana)

# Distances calculation
first_vecino  <- nndist(puntos, k = 1)
second_vecino <- nndist(puntos, k = 2)
third_vecino  <- nndist(puntos, k = 3)
fourth_vecino <- nndist(puntos, k = 4)

# Add distances to db original 
households_distances <- households |>
  mutate(
    first_vecino  = first_vecino,
    second_vecino = second_vecino,
    third_vecino  = third_vecino,
    fourth_vecino = fourth_vecino
    ) |> 
  relocate(geom,.after = fourth_vecino)

write_sf(households_distances,'data/mid/households_distances.gpkg')  