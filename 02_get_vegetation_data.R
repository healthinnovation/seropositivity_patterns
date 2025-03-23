library(rgee)
library(sf)
library(tidyverse)
library(mapview)
library(formattable)
ee_Initialize(quiet = TRUE)

# 1. Reading spatial data -------------------------------------------------
districts <- st_read(
  dsn = 'sources/rawdata/geometry.gpkg',
  layer = 'districts') |> 
  st_simplify(preserveTopology = TRUE, dTolerance = 100)

# 2. spatial data to earth engine object ----------------------------------
ee_districts <- districts |> 
  select(codigo) |> 
  sf_as_ee(quiet = TRUE) 

ee_hydro_ana <- hydro_ana |> 
  select(hydroname) |> 
  sf_as_ee(quiet = TRUE)

ee_hydro_06  <- hydro_06  |> 
  select(hydroname) |> 
  sf_as_ee(quiet = TRUE)

ee_hydro_07  <- hydro_07  |> 
  select(hydroname) |> 
  sf_as_ee(quiet = TRUE)

## 2.1 Parameters
start_date <- 2009
end_date <- 2022

## 2.2 Vegetation variables
ee_leaf_area_index_high  <- 
  ee$ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR') |> 
  ee$ImageCollection$select('leaf_area_index_high_vegetation')

ee_leaf_area_index_high_min <- 
  ee$ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR') |> 
  ee$ImageCollection$select('leaf_area_index_high_vegetation_min')

ee_leaf_area_index_high_max <- 
  ee$ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR') |> 
  ee$ImageCollection$select('leaf_area_index_high_vegetation_max')

ee_leaf_area_index_low  <- 
  ee$ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR') |> 
  ee$ImageCollection$select('leaf_area_index_low_vegetation')

ee_leaf_area_index_low_min <- 
  ee$ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR') |> 
  ee$ImageCollection$select('leaf_area_index_low_vegetation_min')

ee_leaf_area_index_low_max <- 
  ee$ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR') |> 
  ee$ImageCollection$select('leaf_area_index_low_vegetation_max')

ee_evi <- ee$ImageCollection('MODIS/061/MOD13A3') |> 
  ee$ImageCollection$select('EVI')

ee_ndvi <- ee$ImageCollection('MODIS/061/MOD13A3') |> 
  ee$ImageCollection$select('NDVI')

## 2.3 Vegetation variables in R functions
ee_r_leaf_area_index_high <- function(){
  ee_reducer_ic <- ee_leaf_area_index_high$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()
  return(ee_reducer_ic)
}

ee_r_leaf_area_index_high_max <- function(){
  ee_reducer_ic <- ee_leaf_area_index_high_max$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()
  return(ee_reducer_ic)
}

ee_r_leaf_area_index_high_min <- function(){
  ee_reducer_ic <- ee_leaf_area_index_high_min$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()
  return(ee_reducer_ic)
}

ee_r_leaf_area_index_low <- function(){
  ee_reducer_ic <- ee_leaf_area_index_low$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()
  return(ee_reducer_ic)
}

ee_r_leaf_area_index_low_max <- function(){
  ee_reducer_ic <- ee_leaf_area_index_low_max$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()
  return(ee_reducer_ic)
}

ee_r_leaf_area_index_low_min <- function(){
  ee_reducer_ic <- ee_leaf_area_index_low_min$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()
  return(ee_reducer_ic)
}

ee_r_evi <- function(){
  ee_reducer_ic <- ee_evi$filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()$
    multiply(0.0001)
  return(ee_reducer_ic)
}

ee_r_ndvi <- function(){
  ee_reducer_ic <- ee_ndvi$filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()$
    multiply(0.0001)
  return(ee_reducer_ic)
}

# 3. Extracting vegetation variables in the Districts ---------------------

## 3.1 Leaf area index high vegetation
ee_r_leaf_area_index_high() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_leaf_area_index_high

names(districts_leaf_area_index_high) <- sub(
  sprintf("(.{%d})", 7),
  "\\1_",
  gsub('_','',names(districts_leaf_area_index_high)))

## 3.2 Leaf area index high vegetation min
ee_r_leaf_area_index_high_min() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_leaf_area_index_high_min

names(districts_leaf_area_index_high_min) <-  sub(
  sprintf("(.{%d})", 7),
  "\\1_",
  gsub('_','',names(districts_leaf_area_index_high_min)))

## 3.3 Leaf area index high vegetation max
ee_r_leaf_area_index_high_max() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_leaf_area_index_high_max

names(districts_leaf_area_index_high_max) <-  sub(
  sprintf("(.{%d})", 7),
  "\\1_",
  gsub('_','',names(districts_leaf_area_index_high_max)))

## 3.4 Leaf area index low vegetation
ee_r_leaf_area_index_low() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_leaf_area_index_low

names(districts_leaf_area_index_low) <-  sub(
  sprintf("(.{%d})", 7),
  "\\1_",
  gsub('_','',names(districts_leaf_area_index_low)))

## 3.5 Leaf area index low vegetation min
ee_r_leaf_area_index_low_min() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_leaf_area_index_low_min

names(districts_leaf_area_index_low_min) <-  sub(
  sprintf("(.{%d})", 7),
  "\\1_",
  gsub('_','',names(districts_leaf_area_index_low_min)))

## 3.6 Leaf area index low vegetation max
ee_r_leaf_area_index_low_max() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_leaf_area_index_low_max

names(districts_leaf_area_index_low_max) <-  sub(
  sprintf("(.{%d})", 7),
  "\\1_",
  gsub('_','',names(districts_leaf_area_index_low_max)))

## 3.7 EVI
ee_r_evi() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 1000,
    quiet = TRUE) -> districts_evi

names(districts_evi) <-  sub(
  sprintf("(.{%d})", 9),
  "\\1_",
  gsub('_','',names(districts_evi)))

## 3.8 NDVI
ee_r_ndvi() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 1000,
    quiet = TRUE) -> districts_ndvi

names(districts_ndvi) <-  sub(
  sprintf("(.{%d})", 9),
  "\\1_",
  gsub('_','',names(districts_ndvi)))

## 3.9 Dataset final
districts |> 
  left_join(
    y = districts_leaf_area_index_high,
    by = 'codigo') |> 
  left_join(
    y = districts_leaf_area_index_high_min,
    by = 'codigo') |>
  left_join(
    y = districts_leaf_area_index_high_max,
    by = 'codigo') |> 
  left_join(
    y = districts_leaf_area_index_low,
    by = 'codigo') |>
  left_join(
    y = districts_leaf_area_index_low_min,
    by = 'codigo') |>
  left_join(
    y = districts_leaf_area_index_low_max,
    by = 'codigo') |>
  left_join(
    y = districts_evi,
    by = 'codigo') |>
  left_join(
    y = districts_ndvi,
    by = 'codigo') |> 
  st_drop_geometry() |> 
  pivot_longer(
    cols = X200901_leafareaindexhighvegetation:X20221201_NDVI,
    names_to = 'variables',
    values_to = 'valor') |>
  separate(
    col = variables,
    into = c('fecha','variable'),
    sep = '_') |> 
  mutate(
    fecha = str_extract(pattern = '\\d+',fecha),
    fecha = case_when(nchar(fecha) == 6 ~ paste0(fecha,'01'),TRUE ~ fecha),
    fecha = ymd(fecha),
    year = year(fecha), 
    month = month(fecha),
    variable = str_to_lower(variable)) |> 
  relocate(c('variable','valor'),.after = month) |> 
  mutate(variable = str_to_lower(variable)) -> districts_vegetation

if(!dir.exists('output')){dir.create('output')}
if(!dir.exists('output/vegetation')){dir.create('output/vegetation')}
write_csv(districts_vegetation,'output/vegetation/districts_vegetation.csv')