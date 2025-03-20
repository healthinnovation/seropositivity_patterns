library(tidyverse)

# Databases 
ffi_ind <- readRDS("./data/raw/ffi_total.rds") %>% 
  mutate(ffi_h_code =
           str_sub(ffi_is_code, 1, 7)) %>% 
  relocate(ffi_h_code, .after = ffi_is_code)

ffi_hh <- read_csv("./data/raw/household_gps_share_072022.csv") %>% 
  select(ffi_h_code, ffi_gps_lat, ffi_gps_long)

ffi_seropos <- read_csv("./data/raw/FFI Peru_seropositivity_NEG_cutoff.csv") %>% 
  select(Sample, 19:34) %>% 
  mutate(Sample = str_pad(Sample,
                 width = 9,
                 side = "left",
                 pad = "0")) %>% 
  rename(ffi_is_code = Sample)

# Household codes + individual codes + seropositivity

ffi_hh_ind <- ffi_ind %>% 
  full_join(ffi_hh, by = c("ffi_h_code"))

ffi_hh_ind_seropos <-  ffi_hh_ind %>% 
  inner_join(ffi_seropos, by = c("ffi_is_code"))

#write_csv(ffi_hh_ind_seropos, "./data/mid/ffi_hh_ind_seropos.csv")
