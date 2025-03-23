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

# 2. Spatial data to earth engine objects ---------------------------------
ee_districts <- districts |> 
  select(codigo) |> 
  sf_as_ee(quiet = TRUE) 

## 2.1 Parameters
start_date <- 2009
end_date <- 2022

## 2.2 Water variables
ee_aet  <- ee$ImageCollection$Dataset$IDAHO_EPSCOR_TERRACLIMATE |> 
  ee$ImageCollection$select('aet')

ee_def <- ee$ImageCollection$Dataset$IDAHO_EPSCOR_TERRACLIMATE |> 
  ee$ImageCollection$select('def')

ee_pdsi <- ee$ImageCollection$Dataset$IDAHO_EPSCOR_TERRACLIMATE |> 
  ee$ImageCollection$select('pdsi')

ee_pet <- ee$ImageCollection$Dataset$IDAHO_EPSCOR_TERRACLIMATE |> 
  ee$ImageCollection$select('pet')

ee_pr <- ee$ImageCollection$Dataset$IDAHO_EPSCOR_TERRACLIMATE |> 
  ee$ImageCollection$select('pr')

ee_ro <- ee$ImageCollection$Dataset$IDAHO_EPSCOR_TERRACLIMATE |> 
  ee$ImageCollection$select('ro')

ee_soil <- ee$ImageCollection$Dataset$IDAHO_EPSCOR_TERRACLIMATE |> 
  ee$ImageCollection$select('soil')

## 2.3 Water functions 
ee_r_aet <- function(){
  ee_reducer_ic <- ee_aet$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()$
    multiply(0.1)
  return(ee_reducer_ic)
}

ee_r_def <- function(){
  ee_reducer_ic <- ee_def$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()$
    multiply(0.1)
  return(ee_reducer_ic)
}

ee_r_pdsi <- function(){
  ee_reducer_ic <- ee_pdsi$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()$
    multiply(0.01)
  return(ee_reducer_ic)
}

ee_r_pet <- function(){
  ee_reducer_ic <- ee_pet$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()$
    multiply(0.1)
  return(ee_reducer_ic)
}

ee_r_pr <- function(){
  ee_reducer_ic <- ee_pr$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()
  return(ee_reducer_ic)
}

ee_r_ro <- function(){
  ee_reducer_ic <- ee_ro$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()
  return(ee_reducer_ic)
}

ee_r_soil <- function(){
  ee_reducer_ic <- ee_soil$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()$
    multiply(0.1)
  return(ee_reducer_ic)
}

# 3. Extracting water variables in the districts --------------------------

## 3.1 Actual evapotranspiration
ee_r_aet() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_aet

## 3.2 Climate water deficit
ee_r_def() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_def

## 3.3 Palmer Drought Severity Index
ee_r_pdsi() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_pdsi

## 3.4 Reference evapotranspiration
ee_r_pet() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_pet

## 3.5 Precipitation accumulation
ee_r_pr() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_pr

## 3.6 Runoff
ee_r_ro() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_ro

## 3.7 Soil Mousture
ee_r_soil() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_soil


## 3.8 Dataset final
districts |> 
  left_join(
    y = districts_aet,
    by = 'codigo') |> 
  left_join(
    y = districts_def,
    by = 'codigo') |>
  left_join(
    y = districts_pdsi,
    by = 'codigo') |> 
  left_join(
    y = districts_pet,
    by = 'codigo') |>
  left_join(
    y = districts_pr,
    by = 'codigo') |>
  left_join(
    y = districts_ro,
    by = 'codigo') |>
  left_join(
    y = districts_soil,
    by = 'codigo') |> 
  st_drop_geometry() |> 
  pivot_longer(
    cols = X200901_aet:X202212_soil,
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
  relocate(c('variable','valor'),.after = month) -> districts_water

if(!dir.exists('output')){dir.create('output')}
if(!dir.exists('output/water')){dir.create('output/water')}
write_csv(districts_water,'output/water/districts_water.csv')