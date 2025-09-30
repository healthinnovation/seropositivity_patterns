library(tidyverse)
library(sf)
library(lme4)

# Reading datasets --------------------------------------------------------


hotspots <- st_read("./data/mid/gi_star_data.gpkg") %>% 
  st_drop_geometry() %>% 
  select(ffi_h_code, cluster) %>% 
  mutate(
    hot_vs_cold = case_when(
      grepl("Hot spot", cluster) ~ 1L, #Hotspots
      grepl("Cold spot", cluster) ~ 0L, #Coldspots
      TRUE ~ NA_integer_
    ),
    hot_vs_other = case_when(
      grepl("Hot spot", cluster) ~ 1L, #Hotspots
      TRUE ~ 0L #Other categories including not significant
    )
  )

individuals <- read_csv("./data/mid/individuals_multilevel_df.csv") %>% 
  select(-economic_activities) %>% 
  mutate(
    gender = factor(gender,
                    levels = c("Female", "Male")),
    age_oms = factor(age_oms,
                     levels = c("Niño y adolescente", "Adulto joven",
                                "Adulto medio", "Adulto mayor")),
    time_resi_cat2 = factor(time_resi_cat2,
                            levels = c("< 5 años", "5-20 años",
                                       "> 20 años")),
    education_summary = factor(education_summary,
                               levels = c("Sin escolaridad", "Primaria",
                                          "Secundaria",
                                          "Educación superior")),
    sympt_constitucional = factor(sympt_constitucional,
                               levels = c("No", "Sí")),
    sympt_osteomuscular = factor(sympt_osteomuscular,
                                  levels = c("No", "Sí")),
    sympt_gastrointestinal = factor(sympt_gastrointestinal,
                                 levels = c("No", "Sí")),
    sympt_neurologico = factor(sympt_neurologico,
                                    levels = c("No", "Sí")),
    distance_to_rh_hf_sampl = factor(distance_to_rh_hf_sampl,
                                     levels = c("Proximate", "Moderate",
                                                "Distant",
                                                "Extra Distant")),
    EconomicActivities_new = factor(EconomicActivities_new,
                                     levels = c("None", "Student",
                                                "Housewife",
                                                "Farmer", "Fisherman",
                                                "Services/Commerce", "Other")),
    ffi_h_code = factor(ffi_h_code),
    ffi_is_community = factor(ffi_is_community),
    ffi_is_health_facility_name = factor(ffi_is_health_facility_name)
  )


full_multilevel <- individuals %>% 
  full_join(hotspots, by = "ffi_h_code")


df_hot_cold <- full_multilevel %>% 
  drop_na(hot_vs_cold) %>% 
  select(-c(hot_vs_other, cluster))

df_hot_oth <- full_multilevel %>% 
  drop_na(hot_vs_other) %>% 
  select(-c(hot_vs_cold, cluster))

# Multilevel analysis -----------------------------------------------------

hot_cold_empty <- glmer( #Modelo nulo (aca solo esta corriendo la variable de descenlase)
  hot_vs_cold ~ 1 + (1 | ffi_is_health_facility_name/ffi_is_community/ffi_h_code),
  data = df_hot_cold,
  family = poisson(link = "log") #aqui cambiar familia a logistica o a binomial negativa
)

performance::icc(hot_cold_empty) # No estima
sjPlot::tab_model(hot_cold_empty)
broom::tidy(hot_cold_empty)



hot_cold_model <- glmer(
  hot_vs_cold ~ gender + age_oms + time_resi_cat2 + 
  education_summary + sympt_constitucional + sympt_osteomuscular +
  sympt_gastrointestinal + sympt_neurologico +
  distance_to_rh_hf_sampl + EconomicActivities_new +
  (1 | ffi_is_health_facility_name/ffi_is_community),
  data = df_hot_cold,
  family = poisson(link = "log")
)

jtools::summ(hot_cold_model, robust = "HC2", confint = T, digits = 3, exp = T,
                  conf.level = 0.95)


# Hotspots vs Others

hot_oth_empty <- glmer(
  hot_vs_other ~ 1 + (1 | ffi_is_health_facility_name/ffi_is_community/ffi_h_code),
  data = df_hot_oth,
  family = poisson(link = "log")
)

performance::icc(hot_oth_empty) # No estima

hot_oth_model <- glmer(
  hot_vs_other ~ gender + age_oms + time_resi_cat2 + 
    education_summary + sympt_constitucional + sympt_osteomuscular +
    sympt_gastrointestinal + sympt_neurologico +
    distance_to_rh_hf_sampl + EconomicActivities_new +
    (1 | ffi_is_health_facility_name/ffi_is_community/ffi_h_code),
  data = df_hot_oth
)

jtools::summ(hot_cold_model, robust = "HC2", confint = T, digits = 3, exp = T,
             conf.level = 0.95)

# NUEVOS MODELOS----------------------------------------------------------------

# MODELO 1: LOGISTICO
##### MODELO CON HOTS AND COLDS
library(sjPlot)
library(broom.mixed)

hot_cold_empty <- glmer(
  hot_vs_cold ~ 1 + (1 | ffi_is_health_facility_name/ffi_is_community/ffi_h_code),
  data = df_hot_cold,
  family = binomial(link = "logit")
)


performance::icc(hot_cold_empty) # No estima
sjPlot::tab_model(hot_cold_empty)
broom::tidy(hot_cold_empty)

hot_cold_model <- glmer(
  hot_vs_cold ~ gender + age_oms + time_resi_cat2 + 
    education_summary + sympt_constitucional + sympt_osteomuscular +
    sympt_gastrointestinal + sympt_neurologico +
    distance_to_rh_hf_sampl + EconomicActivities_new +
    (1 | ffi_is_health_facility_name/ffi_is_community),
  data = df_hot_cold,
  family = binomial(link = "logit")
)

jtools::summ(hot_cold_model, robust = "HC2", confint = T, digits = 3, exp = T,
             conf.level = 0.95)


# Tabla HTML completa con OR e IC95%
sjPlot::tab_model(
  hot_cold_model,
  transform = "exp",   # convierte log-odds a OR
  show.ci = TRUE,
  show.se = TRUE,
  show.re.var = TRUE,  # incluye varianzas de efectos aleatorios
  show.icc = TRUE,     # muestra ICC si es estimable
  digits = 3
)

# Alternativa: tabla "tidy" en data.frame con OR e IC95%
res_hot_cold <- broom.mixed::tidy(
  hot_cold_model,
  effects = "fixed",
  conf.int = TRUE,
  exponentiate = TRUE
)
print(res_hot_cold)

##### MODELO CON HOTS AND OTHERS
hot_oth_empty <- glmer(
  hot_vs_other ~ 1 + (1 | ffi_is_health_facility_name/ffi_is_community/ffi_h_code),
  data = df_hot_oth,
  family = binomial(link = "logit")
)

performance::icc(hot_oth_empty) # No estima

hot_oth_model <- glmer(
  hot_vs_other ~ gender + age_oms + time_resi_cat2 + 
    education_summary + sympt_constitucional + sympt_osteomuscular +
    sympt_gastrointestinal + sympt_neurologico +
    distance_to_rh_hf_sampl + EconomicActivities_new +
    (1 | ffi_is_health_facility_name/ffi_is_community/ffi_h_code),
  data = df_hot_oth,
  family = binomial(link = "logit")
)

jtools::summ(hot_oth_model, robust = "HC2", confint = T, digits = 3, exp = T,
             conf.level = 0.95)

sjPlot::tab_model(
  hot_oth_model,
  transform = "exp",   # convierte log-odds a Odds Ratios
  show.ci = TRUE,
  show.se = TRUE,
  show.re.var = TRUE,  # incluye varianzas de efectos aleatorios
  show.icc = TRUE,     # muestra ICC si es estimable
  digits = 3
)

# ============================
# TABLA "TIDY" (data frame con OR e IC95%)
# ============================
res_hot_oth <- broom.mixed::tidy(
  hot_oth_model,
  effects = "fixed",      # solo efectos fijos
  conf.int = TRUE,        # intervalos de confianza
  exponentiate = TRUE     # convierte a OR
)

print(res_hot_oth)
####################### MODELO CON BINOMIAL NEGATIVO

library(lme4)
library(performance)
library(sjPlot)
library(broom.mixed)

# ============================
# MODELO NULO NB
# ============================
hot_cold_empty_nb <- glmer.nb(
  hot_vs_cold ~ 1 + (1 | ffi_is_health_facility_name/ffi_is_community/ffi_h_code),
  data = df_hot_cold,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)
performance::icc(hot_cold_empty_nb) # No estima

hot_cold_empty_nb <- glmer.nb(
  hot_vs_other ~ gender + age_oms + time_resi_cat2 + 
    education_summary + sympt_constitucional + sympt_osteomuscular +
    sympt_gastrointestinal + sympt_neurologico +
    distance_to_rh_hf_sampl + EconomicActivities_new +
    (1 | ffi_is_health_facility_name/ffi_is_community/ffi_h_code),
  data = df_hot_oth,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

jtools::summ(hot_cold_empty_nb, robust = "HC2", confint = TRUE, digits = 3, exp = TRUE,
             conf.level = 0.95)

sjPlot::tab_model(
  hot_cold_empty_nb,
  transform = "exp",
  show.ci = TRUE,
  show.se = TRUE,
  show.re.var = TRUE,
  show.icc = TRUE,
  digits = 3
)

res_hot_oth_nb <- broom.mixed::tidy(
  hot_cold_empty_nb,
  effects = "fixed",
  conf.int = TRUE,
  exponentiate = TRUE
)

print(res_hot_oth_nb)


# ===== HOT vs OTHER — Binomial Negativo =====

# Modelo nulo NB
hot_oth_empty_nb <- glmer.nb(
  hot_vs_other ~ 1 + (1 | ffi_is_health_facility_name/ffi_is_community/ffi_h_code),
  data   = df_hot_oth,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

performance::icc(hot_oth_empty_nb) # Puede no estimar; si falla, usar r2_nakagawa()

# Modelo con predictores NB
hot_oth_model <- glmer.nb(
  hot_vs_other ~ gender + age_oms + time_resi_cat2 + 
    education_summary + sympt_constitucional + sympt_osteomuscular +
    sympt_gastrointestinal + sympt_neurologico +
    distance_to_rh_hf_sampl + EconomicActivities_new +
    (1 | ffi_is_health_facility_name/ffi_is_community/ffi_h_code),
  data   = df_hot_oth,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

# Resumen (razones de tasas exp(beta))
jtools::summ(hot_oth_model, robust = TRUE, confint = TRUE, digits = 3, exp = TRUE,
             conf.level = 0.95)

sjPlot::tab_model(
  hot_oth_model,
  transform = "exp",
  show.ci = TRUE,
  show.se = TRUE,
  show.re.var = TRUE,
  show.icc = TRUE,
  digits = 3
)

res_hot_oth_nb <- broom.mixed::tidy(
  hot_oth_model,
  effects = "fixed",
  conf.int = TRUE,
  exponentiate = TRUE
)
print(res_hot_oth_nb)

