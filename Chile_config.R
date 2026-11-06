suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(parallel)
})

# Configuración común a todos los scripts del análisis de Chile. La cargan
# Chile_comun.R (preparación completa) y los scripts que solo leen salidas ya
# calculadas. Es el único sitio donde se cambian rutas, escala y umbral.

# 0 Configuración ----

DIR_BASE        <- "C:/A_TRABAJO/CHILE"
DIR_CHELSA      <- "C:/A_TRABAJO/DATA/CHELSA/BIOCLIM_1981_2010"
PATH_STUDY_AREA <- file.path(DIR_BASE, "DATA/Limites_admin/chile_limites_10km_diss.shp")
PATH_REGIONES   <- file.path(DIR_BASE, "DATA/Limites_admin/chile_regiones_10km.shp")
PATH_MACROZONAS <- file.path(DIR_BASE, "DATA/Limites_admin/chile_macrozona_10km.shp")
COL_REGION      <- "REGION"
COL_MACROZONA   <- "MACROZONA"
PATH_KOPPEN     <- file.path(DIR_BASE, "DATA/koppen/layer_koppen/layer_koppen_20220309023154.shp")
COL_KOPPEN      <- NULL       # campo con el tipo de clima; NULL lo detecta
PATH_AP_RAW     <- file.path(DIR_BASE, "DATA/Areas_Protegidas/APs_chile_2026.shp")
PATH_AP         <- file.path(DIR_BASE, "DATA/Areas_Protegidas/APs_chile_2026_5km2.shp")

COL_ID      <- "SITECODE"
TH          <- 0.99       # umbral de las envolventes
VIF_TH      <- 5          # solo diagnóstico, al construir el clima
MIN_AREA_M2 <- 5e6        # 5 km2, tras recortar con el área de estudio
UNIT_KM     <- 10         # lado de la malla de referencia
RES_EA      <- 1000       # resolución de la rejilla equivalente final (m)
NUM_CORES   <- 8
RCRI_LOG    <- log10

# LAEA centrada en Chile continental, elipsoide GRS80 (SIRGAS).
CRS_EA <- "+proj=laea +lat_0=-37 +lon_0=-71 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"

P        <- round(TH * 100)
DIR_RES  <- file.path(DIR_BASE, "RESULTADOS")
# Escenario. 0 es la red actual; un valor mayor añade a cada AP un anillo de
# BUFFER_KM km que excluye el terreno de cualquier AP (zona de amortiguación).
BUFFER_KM <- 0
dir_escenario <- function(buffer_km)
  file.path(DIR_RES, sprintf("APCRS_P%d_%dkm%s", P, UNIT_KM, if (buffer_km > 0) sprintf("_B%gkm", buffer_km) else ""))
DIR_OUT  <- dir_escenario(BUFFER_KM)
TAG_ESC  <- if (BUFFER_KM > 0) sprintf("_B%gkm", BUFFER_KM) else ""   # sufijo de escenario en los nombres de salida
DIR_TMP  <- file.path(DIR_OUT, "tmp")
CLIM_VIF <- file.path(DIR_RES, "clima", "present_climate_VIF_WGS84.tif")
TAG_BA   <- sprintf("BCRS_area_%dkm_P%d", UNIT_KM, P)
TAG_AA   <- sprintf("APCRS_area_%dkm_P%d", UNIT_KM, P)
TAG_RA   <- sprintf("RCRI_area_P%d", P)
for (d in c(DIR_TMP, dirname(CLIM_VIF))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
terraOptions(tempdir = DIR_TMP)
sf_use_s2(FALSE)
s_m <- UNIT_KM * 1000
# Capa de AP del escenario (la red actual o las AP con su anillo).
PATH_AP_ESC <- if (BUFFER_KM > 0) file.path(DIR_OUT, "APs_escenario.gpkg") else PATH_AP
message("Escala ", UNIT_KM, " km | P", P, " | anillo ", BUFFER_KM, " km | salida ", DIR_OUT)

RES_RUN <- list(UNIT_KM = UNIT_KM, TH = TH, BUFFER_KM = BUFFER_KM)
anota <- function(...) { x <- list(...); RES_RUN[names(x)] <<- x; invisible(NULL) }
