library(rgee)
library(sf)
library(tidyverse)
library(lubridate)
ee_Initialize(quiet = TRUE)
sf_use_s2(use_s2 = F)

# 1. Reading spatial data -------------------------------------------------
districts <- read_csv(file = 'data/mid/ffi_hh_ind_seropos.csv') |> 
  st_transform(crs = 32718) |> 
  st_simplify(preserveTopology = TRUE, dTolerance = 500) |> 
  st_transform(crs = 4326)

# 2. Spatial data to earth engine objects ---------------------------------
ee_districts <- districts|> 
  dplyr::select(codigo)  |> 
  sf_as_ee(quiet = TRUE) 

ee_hydro_ana <- hydro_ana |> 
  dplyr::select(hydroname)|> 
  sf_as_ee(quiet = TRUE)

ee_hydro_06  <- hydro_06  |> 
  dplyr::select(hydroname)|>
  sf_as_ee(quiet = TRUE)

ee_hydro_07  <- hydro_07  |> 
  dplyr::select(hydroname)|> 
  sf_as_ee(quiet = TRUE)

## 2.1 Parameters
start_date <- 2009
end_date   <- 2022

ee_tmin <- ee$ImageCollection$Dataset$IDAHO_EPSCOR_TERRACLIMATE |> 
  ee$ImageCollection$select('tmmn')

ee_tmax <- ee$ImageCollection$Dataset$IDAHO_EPSCOR_TERRACLIMATE |> 
  ee$ImageCollection$select('tmmx')

ee_wsp <- ee$ImageCollection$Dataset$IDAHO_EPSCOR_TERRACLIMATE |> 
  ee$ImageCollection$select('vs')

ee_temp_2m <- ee$ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR') |> 
  ee$ImageCollection$select('temperature_2m')

ee_temp_soil_7cm_deep <- ee$ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR') |> 
  ee$ImageCollection$select('soil_temperature_level_1')

ee_etp_max <- ee$ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR') |> 
  ee$ImageCollection$select('potential_evaporation_max')

ee_etp_min <- ee$ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR') |> 
  ee$ImageCollection$select('potential_evaporation_min')

ee_water_rainfall <- ee$ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR') |> 
  ee$ImageCollection$select('runoff_sum')

ee_accu_liquid <- ee$ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR') |> 
  ee$ImageCollection$select('total_precipitation_sum')

ee_day_lst <- ee$ImageCollection('MODIS/061/MYD21C3')$select('LST_Day')

## 2.2 Climate variables in R functions
ee_r_tmin <- function(){
  ee_reducer_ic <- ee_tmin$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()$
    multiply(0.1)
  return(ee_reducer_ic)
}

ee_r_tmax <- function(){
  ee_reducer_ic <- ee_tmax$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()$
    multiply(0.1)
  return(ee_reducer_ic)
}

ee_r_wsp <- function(){
  ee_reducer_ic <- ee_wsp$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()$
    multiply(0.01)
  return(ee_reducer_ic)
}

ee_r_temp_2m <- function(){
  ee_reducer_ic <- ee_temp_2m$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands() |> 
    ee$Image$subtract(273.15)
  return(ee_reducer_ic)
}

ee_r_temp_soil_7cm_deep <- function(){
  ee_reducer_ic <- ee_temp_soil_7cm_deep$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands() |> 
    ee$Image$subtract(273.15)
  return(ee_reducer_ic)
}

ee_r_etp_max <-  function(){
  ee_reducer_ic <- ee_etp_max$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()
  return(ee_reducer_ic)
}


ee_r_etp_min <-  function(){
  ee_reducer_ic <- ee_etp_min$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()
  return(ee_reducer_ic)
}


ee_r_water_rainfall <- function(){
  ee_reducer_ic <- ee_water_rainfall$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()
  return(ee_reducer_ic)
}

ee_r_accu_liquid <- function(){
  ee_reducer_ic <- ee_accu_liquid$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands()
  return(ee_reducer_ic)
}

ee_r_day_lst <- function(){
  ee_reducer_ic <- ee_day_lst$
    filter(ee$Filter$calendarRange(start_date,end_date,'year'))$
    toBands() |> 
    ee$Image$subtract(273.15)
  return(ee_reducer_ic)
}

# 3. Extracting climate variables in the districts ------------------------

## 3.1 Temperature minimum
ee_r_tmin() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$min(),
    sf = FALSE,
    scale = 5000,
    quiet = TRUE) -> districts_tmin

## 3.2 Temperature maximum
ee_r_tmax() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$max(),
    sf = FALSE,
    scale = 5000,
    quiet = TRUE) -> districts_tmax


## 3.3 Wind speed
ee_r_wsp() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 5000,
    quiet = TRUE) -> districts_wsp


## 3.4 Temperature of air at 2m
ee_r_temp_2m() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_temp_2m

names(districts_temp_2m) <- str_c(
  str_sub(names(districts_temp_2m), 1, 19),
  str_sub(names(districts_temp_2m), 21))

## 3.5 Temperature of the soil in layer (0-7cm)
ee_r_temp_soil_7cm_deep() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_temp_soil_7cm_deep

names(districts_temp_soil_7cm_deep) <- str_replace(
  names(districts_temp_soil_7cm_deep),
  '_temperature_level_1$',
  'temperature')

## 3.6 Potential evaporation max
ee_r_etp_max() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$max(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_etp_max

names(districts_etp_max) <- str_replace(
  names(districts_etp_max),
  "potential_evaporation_max$",
  "potentialevaporation.max")

## 3.7 Potential evaporation min
ee_r_etp_min() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$min(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_etp_min

names(districts_etp_min) <- str_replace(
  names(districts_etp_min),
  "potential_evaporation_min$",
  "potentialevaporation.min")

## 3.8 Some water from rainfall
ee_r_water_rainfall() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_water_rainfall

names(districts_water_rainfall) <- str_replace(
  names(districts_water_rainfall),
  "runoff_sum$",
  "runoff")

## 3.9 Accumulated liquid and frozen water, including rain and snow, that falls to the Earth's surface
ee_r_accu_liquid() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$sum(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_accu_liquid

names(districts_accu_liquid) <- str_replace(
  names(districts_accu_liquid),
  "total_precipitation_sum$",
  "totalprecipitation")

## 3.10 Daytime Land Surface Temperature
ee_r_day_lst() |> 
  ee_extract(
    y = ee_districts,
    fun = ee$Reducer$mean(),
    sf = FALSE,
    scale = 11000,
    quiet = TRUE) -> districts_day_lst

names(districts_day_lst) <- sub(
  "^(X\\d{4})_(\\d{2})_(\\d{2})_(LST)_(Day)$",
  "\\1\\2\\3_\\4\\5",
  names(districts_day_lst))


## 3.11 Dataset of climate factor for districts
districts |> 
  left_join(
    y = districts_tmin,
    by = 'codigo') |> 
  left_join(
    y = districts_tmax,
    by = 'codigo') |>
  left_join(
    y = districts_wsp,
    by = 'codigo') |> 
  left_join(
    y = districts_temp_2m,
    by = 'codigo') |>
  left_join(
    y = districts_temp_soil_7cm_deep,
    by = 'codigo') |>
  left_join(
    y = districts_etp_max,
    by = 'codigo') |>
  left_join(
    y = districts_etp_min,
    by = 'codigo') |>
  left_join(
    y = districts_water_rainfall,
    by = 'codigo') |>
  left_join(
    y = districts_accu_liquid,
    by = 'codigo') |>
  left_join(
    y = districts_day_lst,
    by = 'codigo') |> 
  st_drop_geometry() |> 
  pivot_longer(
    cols = X200901_tmmn:X20221201_LSTDay,
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
    month = month(fecha)) |> 
  relocate(c('variable','valor'),.after = month) |> 
  mutate(variable = str_to_lower(variable)) -> districts_climate

if(!dir.exists('output')){dir.create('output')}
if(!dir.exists('output/climate')){dir.create('output/climate')}
write_csv(districts_climate,'output/climate/districts_climate.csv')