# Efecto del cambio climático sobre varias AP, como ejemplo para la
# presentación. Reserva Nacional Coyhaique (60RN) y Parque Nacional Cerro
# Castillo (18PN) en Aysén, y Parques Nacionales La Campana (29PN) y Río
# Clarillo (48PN) en la zona central. Los nombres de las salidas conservan el
# sufijo Aysen para no romper rutas ya usadas.
#
# Dos pasos.
#   1 mh_rep_ch() de ClimaRep con el clima presente del análisis de Chile y un
#     clima futuro preparado aquí con el mismo tratamiento. Clasifica cada celda
#     de Chile como no análoga, estable, perdida o nueva respecto al clima de
#     cada AP.
#   2 Curva de aislamiento a partir de las salidas Change, Mh_Raw_Pre y
#     Mh_Raw_Fut, con los índices dT y L, la bandera de ausencia de sustituto y
#     los descriptores (A, T_pre, S_int, S_fut, D_min, D_eq, H, W), siguiendo
#     el método de AISLAMIENTO_TEST_15_9.R.
#
# Diferencias con AISLAMIENTO_TEST_15_9.R, y por qué.
#   a) Superficie en km2 y no en celdas. mh_rep_ch exige un CRS con código EPSG
#      y el clima del análisis está en WGS84, así que las salidas quedan en
#      celdas de 30 segundos, cuya superficie cambia con la latitud (unos 0.8
#      km2 a 20 S frente a 0.55 km2 a 50 S). Cada celda pesa por su superficie
#      real (cellSize), sin remuestrear. En una rejilla de igual área es lo
#      mismo que contar celdas. El número de celdas se conserva para el control
#      de tamaño mínimo (A_MINIMO).
#   b) El radio del envolvente se reconstruye con las celdas por centro. La
#      versión de mh_rep_ch de este proyecto usa mask(touches = FALSE) para
#      definir las celdas de cada polígono (y touches = TRUE solo si ningún
#      centro cae dentro). La nota 8c del script de aislamiento, que usaba la
#      muestra por contacto, corresponde a una versión anterior del paquete. Ya
#      no hay inflación de borde que medir y el control del radio debe cuadrar
#      exactamente.
#   c) Pocas AP, así que no hay estratos de H ni correlaciones entre AP.
#   d) Rótulos en español, porque el destino es la presentación a gestores.
#
# Umbral. [MML, 6 nov 2026] P99, por coherencia con el análisis presente de
# Chile. Se aparta del cuantil canónico 0.95 del método de aislamiento ibérico,
# así que estas cifras no son comparables con las de los parques nacionales
# ibéricos. La banda de sensibilidad de la figura se lee a 0.95 y 0.995.


# 0 Configuración ----

PATH_CONFIG    <- "Chile_config.R"   # configuración común del análisis de Chile
PATH_MH_REP_CH <- "mh_rep_ch.R"      # versión de desarrollo con cov_type; NULL usa ClimaRep instalado
source(PATH_CONFIG, encoding = "UTF-8")
suppressPackageStartupMessages({
  library(ggplot2)     # mh_rep_ch usa element_text() sin prefijo cuando se carga con source()
  library(tidyterra)
})

AP_EJEMPLO <- c("60RN", "18PN", "29PN", "48PN")
NOMBRES_AP <- c("60RN" = "Reserva Nacional Coyhaique",
                "18PN" = "Parque Nacional Cerro Castillo",
                "29PN" = "Parque Nacional La Campana",
                "48PN" = "Parque Nacional Río Clarillo")

# Clima futuro. [MML] CHELSA V2.1, GFDL-ESM4, SSP5-8.5, 2041-2070. Archivos del
# tipo CHELSA_bio1_2041-2070_gfdl-esm4_ssp585_V.2.1.tif. PERIODO debe escribirse
# como en el nombre del archivo (2041-2070), no como en la carpeta (2040_2070).
# Con varios GCM se usaría la media de cada variable (ensemble).
DIR_CHELSA_FUT <- "C:/A_TRABAJO/DATA/CHELSA/FUTURE/2041_2071/SSP585/IPSL"
GCMS    <-  "ipsl-cm6a-lr" #"gfdl-esm4"
SSP     <- "ssp585"
PERIODO <- "2041-2070"

TH_CH    <- TH                       # 0.99 de Chile_config.R, el mismo P99 del análisis presente
TH_BANDA <- if (TH_CH < 0.99) c(0.90, TH_CH, 0.99) else c(0.95, TH_CH, 0.995)
COV_TYPE <- "present"                # huella presente idéntica a la de mh_rep()
stopifnot(TH_CH == 0.99, length(TH_BANDA) == 3, !is.unsorted(TH_BANDA), TH_BANDA[2] == TH_CH)

ANCHO_BANDA_KM <- 1                  # solo para dibujar
N_VARIABLES    <- 5
A_MINIMO       <- 10 * N_VARIABLES   # celdas, diez por dimensión
MARGEN_MAPA_KM <- 120                # margen alrededor de cada AP en su mapa

MODELO   <- if (length(GCMS) > 1) paste0("ENSEMBLE-", SSP) else paste0(GCMS, "-", SSP)
TAG_FUT  <- paste(MODELO, PERIODO, sep = "_")
DIR_CH   <- file.path(DIR_RES, "CAMBIO_AYSEN", sprintf("%s_P%d", TAG_FUT, round(100 * TH_CH)))
DIR_FIG  <- file.path(DIR_CH, "figuras")
CLIM_FUT <- file.path(DIR_RES, "clima", sprintf("future_climate_%s_WGS84.tif", TAG_FUT))
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)

ETIQ_ESC <- sprintf("%s, %s, %s", if (length(GCMS) > 1) sprintf("media de %d GCM", length(GCMS)) else GCMS,
                    toupper(SSP), PERIODO)
message("Escenario futuro: ", ETIQ_ESC, " | umbral ", TH_CH, " | salida ", DIR_CH)
t0 <- Sys.time()


# 1 Clima presente y futuro ----
# El futuro se lee y se recorta exactamente igual que el presente en
# Chile_comun.R (mismas variables sin transformar, mismo CRS, mismo recorte con
# el área de estudio). Se comprueba después que presente y futuro estén en la
# misma escala, porque una diferencia de escala u offset entre archivos CHELSA
# arruinaría las distancias sin dar ningún error.

VARS <- c("bio1", "bio2", "bio4", "bio12", "bio15")
clim <- rast(CLIM_VIF)
stopifnot(identical(names(clim), VARS), crs(clim, describe = TRUE)$code == "4326")

sa <- read_sf(PATH_STUDY_AREA) |> st_make_valid() |> st_transform(4326) |> st_union() |> st_make_valid()
sa <- st_sf(ID = 1L, geometry = sa)
sa_v <- vect(sa)

var_de <- function(f)
  sub("^bio0+([1-9])", "bio\\1", sub("^.*?(bio[0-9]+).*$", "\\1", tolower(basename(f)), perl = TRUE))

f_todos <- list.files(DIR_CHELSA_FUT, "\\.tif$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
ficheros_gcm <- function(gcm) {
  b <- tolower(basename(f_todos))
  f <- f_todos[grepl(gcm, b, fixed = TRUE) & grepl(SSP, b, fixed = TRUE) &
                 grepl(PERIODO, b, fixed = TRUE) & grepl("bio[0-9]+", b)]
  v <- var_de(f)
  if (anyDuplicated(v[v %in% VARS])) stop("Variables repetidas para ", gcm, ". Revisar ", DIR_CHELSA_FUT)
  i <- match(VARS, v)
  if (anyNA(i)) stop("Faltan variables para ", gcm, ": ", paste(VARS[is.na(i)], collapse = ", "))
  f[i]
}
fuentes  <- unlist(lapply(GCMS, ficheros_gcm))
f_fuente <- sub("\\.tif$", "_fuentes.csv", CLIM_FUT)

# Se reutiliza el futuro ya construido solo si procede de los mismos archivos.
reutiliza <- file.exists(CLIM_FUT) && file.exists(f_fuente) &&
  identical(sort(basename(fuentes)), sort(read.csv(f_fuente)$fuente))
if (!reutiliza) {
  capas <- lapply(GCMS, function(g) {
    r <- rast(ficheros_gcm(g)); names(r) <- VARS; crs(r) <- "EPSG:4326"
    r <- mask(crop(r, sa_v, snap = "out"), sa_v, touches = FALSE)
    if (!compareGeom(r, clim, stopOnError = FALSE))
      stop("La rejilla futura de ", g, " no coincide con la del presente. No se remuestrea en silencio.")
    r
  })
  clim_fut <- Reduce(`+`, capas) / length(capas)
  names(clim_fut) <- VARS
  writeRaster(clim_fut, CLIM_FUT, overwrite = TRUE,
              gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "BIGTIFF=IF_SAFER"))
  write.csv(data.frame(fuente = basename(fuentes)), f_fuente, row.names = FALSE)
  rm(capas)
}
clim_fut <- rast(CLIM_FUT)
stopifnot(identical(names(clim_fut), VARS), compareGeom(clim, clim_fut))

g_pre <- global(clim, c("mean", "sd"), na.rm = TRUE)
g_fut <- global(clim_fut, c("mean", "sd"), na.rm = TRUE)
escala <- data.frame(var = VARS, media_pre = g_pre$mean, media_fut = g_fut$mean,
                     sd_pre = g_pre$sd, sd_fut = g_fut$sd,
                     cambio_en_sd = (g_fut$mean - g_pre$mean) / g_pre$sd,
                     razon_sd = g_fut$sd / g_pre$sd)
message("\nPresente y futuro en el área de estudio")
print(escala, digits = 3)
if (any(abs(escala$cambio_en_sd) > 1.5 | escala$razon_sd < 0.5 | escala$razon_sd > 2))
  stop("Presente y futuro no parecen estar en la misma escala (factor u offset distinto). Revisar los archivos CHELSA.")
n_pre <- global(!is.na(clim[[1]]), "sum")[[1]]
n_fut <- global(!is.na(clim_fut[[1]]), "sum")[[1]]
message(sprintf("Celdas con clima presente %d | futuro %d", n_pre, n_fut))
write.csv(escala, file.path(DIR_CH, "control_escala_presente_futuro.csv"), row.names = FALSE)


# 2 Áreas protegidas y región ----
# Las mismas geometrías que en el análisis presente (agrupadas por SITECODE,
# recortadas con el área de estudio y de más de 5 km2).

aps_all <- read_sf(PATH_AP) |> st_make_valid() |> st_transform(4326)
stopifnot(COL_ID %in% names(aps_all))
aps <- aps_all[aps_all[[COL_ID]] %in% AP_EJEMPLO, ]
if (nrow(aps) != length(AP_EJEMPLO))
  stop("No están todas las AP en ", PATH_AP, ": ", paste(setdiff(AP_EJEMPLO, aps[[COL_ID]]), collapse = ", "))
aps <- aps[match(AP_EJEMPLO, aps[[COL_ID]]), COL_ID]
aps$nombre   <- unname(NOMBRES_AP[aps[[COL_ID]]])
aps$area_km2 <- as.numeric(st_area(st_transform(aps, CRS_EA))) / 1e6

# Región de cada AP, la de mayor superficie compartida (en LAEA). El resumen
# regional se hace en esa región. Las celdas de costa que caen fuera de la capa
# de regiones no cuentan en él.
reg    <- read_sf(PATH_REGIONES) |> st_make_valid() |> st_transform(4326)
reg_ea <- st_transform(reg, CRS_EA); aps_ea <- st_transform(aps, CRS_EA)
ov <- st_intersects(aps_ea, reg_ea)
aps$region <- vapply(seq_len(nrow(aps)), function(i) {
  if (!length(ov[[i]])) stop("La AP ", aps[[COL_ID]][i], " no cae en ninguna región.")
  a <- vapply(ov[[i]], function(j)
    sum(as.numeric(st_area(suppressWarnings(st_intersection(st_geometry(aps_ea)[i], st_geometry(reg_ea)[j]))))),
    numeric(1))
  as.character(reg[[COL_REGION]][ov[[i]][which.max(a)]])
}, character(1))
reg_ap <- lapply(aps$region, function(nm) st_sf(geometry = st_union(st_geometry(reg[reg[[COL_REGION]] == nm, ]))))
rm(reg_ea, aps_ea, ov)
print(st_drop_geometry(aps))


# 3 mh_rep_ch ----
# Se borran antes las salidas esperadas, para que nada proceda de una ejecución
# anterior.

fun_mh <- if (!is.null(PATH_MH_REP_CH)) {
  source(PATH_MH_REP_CH, encoding = "UTF-8", local = TRUE); mh_rep_ch
} else ClimaRep::mh_rep_ch
stopifnot("cov_type" %in% names(formals(fun_mh)))

slug <- function(x) {          # mismo nombre de archivo que construye mh_rep_ch
  x <- iconv(tolower(as.character(x)), to = "ASCII//TRANSLIT")
  x <- gsub("[^a-z0-9_]+", "_", x); x <- gsub("__+", "_", x); gsub("^_|_$", "", x)
}
aps$slug <- slug(aps[[COL_ID]])
f_ch <- file.path(DIR_CH, "Change",     paste0("Th_change_",  MODELO, "_", PERIODO, "_", aps$slug, ".tif"))
f_mp <- file.path(DIR_CH, "Mh_Raw_Pre", paste0("Mh_raw_Pre_", aps$slug, ".tif"))
f_mf <- file.path(DIR_CH, "Mh_Raw_Fut", paste0("Mh_raw_Fut_", MODELO, "_", PERIODO, "_", aps$slug, ".tif"))
unlink(c(f_ch, f_mp, f_mf))

fun_mh(polygon = aps, col_name = COL_ID,
       present_climate_variables = clim, future_climate_variables = clim_fut,
       study_area = sa, th = TH_CH, model = MODELO, year = PERIODO,
       dir_output = DIR_CH, save_raw = TRUE, cov_type = COV_TYPE)
if (!all(file.exists(c(f_ch, f_mp, f_mf)))) stop("mh_rep_ch no ha generado todas las salidas esperadas.")
message("mh_rep_ch listo | ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")


# 4 Curva de aislamiento ----

# Distancia a la que un conjunto de celdas acumula 'objetivo' km2. Exacta, sin
# interpolar dentro de bandas.
dist_area <- function(d, w, objetivo) {
  if (!length(d) || !is.finite(objetivo) || objetivo <= 0) return(NA_real_)
  o <- order(d); i <- which(cumsum(w[o]) >= objetivo * (1 - 1e-9))[1]
  if (is.na(i)) NA_real_ else d[o][i]
}
suma_banda <- function(b, w, NB) {
  s <- numeric(NB)
  if (length(b)) { r <- rowsum(w, b); s[as.integer(rownames(r))] <- r[, 1] }
  s
}
# Superficie por banda de distancia. Solo para dibujar.
matriz_bandas <- function(d, w, sel, ancho) {
  NB <- max(1L, ceiling(max(d) / ancho))
  b  <- pmax(1L, ceiling(d / ancho))
  out <- data.frame(banda = seq_len(NB), d_km = seq_len(NB) * ancho, disponibles = suma_banda(b, w, NB))
  for (nm in names(sel)) out[[nm]] <- suma_banda(b[sel[[nm]]], w[sel[[nm]]], NB)
  out[seq_len(max(which(out$disponibles > 0), 1L)), ]
}

procesar_ap <- function(i) {
  codigo <- aps[[COL_ID]][i]
  r    <- rast(f_ch[i]); r_mp <- rast(f_mp[i]); r_mf <- rast(f_mf[i])
  if (!compareGeom(r, r_mp, stopOnError = FALSE) || !compareGeom(r, r_mf, stopOnError = FALSE))
    stop("Change y Mahalanobis no comparten rejilla en ", codigo)
  ap_v <- vect(aps[i, ])
  m_i  <- rasterize(ap_v, r)                                  # centro de celda
  v_i  <- !is.na(values(m_i, mat = FALSE))
  v_r  <- !is.na(values(rasterize(vect(reg_ap[[i]]), r), mat = FALSE))
  v_c  <- values(r, mat = FALSE)
  v_a  <- values(cellSize(r, unit = "km"), mat = FALSE)
  v_d  <- values(distance(m_i), mat = FALSE) / 1000
  vp   <- values(r_mp, mat = FALSE); vf <- values(r_mf, mat = FALSE)
  if (!any(v_i)) stop("Ninguna celda con centro dentro de ", codigo)

  ok     <- !is.na(v_c)
  an_pre <- ok & v_c %in% c(1, 2)                             # estable + perdida
  an_fut <- ok & v_c %in% c(1, 3)                             # estable + nueva

  # Clasificación canónica. A es la superficie interior análoga a su propio
  # clima y queda fija como listón.
  A      <- sum(v_a[v_i & an_pre]); n_A <- sum(v_i & an_pre)
  if (A == 0) stop("A = 0 en ", codigo)
  int_es <- sum(v_a[v_i & ok & v_c == 1]); int_pe <- sum(v_a[v_i & ok & v_c == 2])
  int_nu <- sum(v_a[v_i & ok & v_c == 3])
  S_int  <- (int_es + int_nu) / A

  ext   <- !v_i & ok & is.finite(v_d)                         # territorio exterior
  a_ext <- sum(v_a[ext])
  a_pre <- sum(v_a[ext & an_pre]); a_fut <- sum(v_a[ext & an_fut])
  T_pre <- a_pre / a_ext; T_fut <- a_fut / a_ext
  S_pre <- a_pre / A;     S_fut <- a_fut / A
  dT    <- if (a_pre > 0) a_fut / a_pre else NA_real_

  D_min    <- dist_area(v_d[ext], v_a[ext], A)
  D_eq     <- dist_area(v_d[ext & an_fut], v_a[ext & an_fut], A)
  D_eq_pre <- dist_area(v_d[ext & an_pre], v_a[ext & an_pre], A)
  bandera  <- as.integer(S_fut < 1)                           # tiene precedencia sobre L
  L <- if (bandera == 0 && is.finite(D_min) && D_min > 0 && is.finite(D_eq)) D_eq / D_min else NA_real_

  # Región de la AP, fuera de la propia AP.
  er <- ext & v_r
  reg_pre <- sum(v_a[er & an_pre]); reg_fut <- sum(v_a[er & an_fut])
  reg_est <- sum(v_a[er & ok & v_c == 1])

  M_b <- matriz_bandas(v_d[ext], v_a[ext],
                       list(stable = (v_c == 1)[ext], lost = (v_c == 2)[ext], novel = (v_c == 3)[ext]),
                       ANCHO_BANDA_KM)
  M_b$ap <- codigo

  # Envolvente. Radios sobre las celdas interiores por centro con presente
  # válido, igual que mh_rep_ch. c_lo y c_hi encierran el radio del paquete.
  fin_p <- is.finite(vp); fin <- fin_p & is.finite(vf)
  x_int <- vp[v_i & fin_p]
  H <- mean(x_int); H_med <- stats::median(x_int)
  cc <- as.numeric(stats::quantile(x_int, probs = TH_BANDA))
  W  <- cc[3] / cc[1]
  c_lo <- max(vp[an_pre & fin_p]); c_hi <- min(vp[ok & !(v_c %in% c(1, 2)) & fin_p])

  # Banda de sensibilidad. Cada curva va normalizada por su propio listón.
  ext_m <- ext & fin
  d_m <- v_d[ext_m]; vf_m <- vf[ext_m]; a_m <- v_a[ext_m]
  A_th <- S_th <- rep(NA_real_, 3); sel <- list()
  for (k in 1:3) {
    A_th[k] <- sum(v_a[v_i & fin & vp <= cc[k]])
    sel[[c("fut_lo", "fut_mid", "fut_hi")[k]]] <- vf_m <= cc[k]
    S_th[k] <- if (A_th[k] > 0) sum(a_m[vf_m <= cc[k]]) / A_th[k] else NA_real_
  }
  M_th <- matriz_bandas(d_m, a_m, sel, ANCHO_BANDA_KM); M_th$ap <- codigo

  fila <- data.frame(
    ap = codigo, nombre = aps$nombre[i], area_km2 = aps$area_km2[i],
    dT = dT, L = L, bandera = bandera,
    A = A, n_A = n_A, n_int = sum(v_i & fin_p), frac_A = n_A / sum(v_i & ok),
    T_pre = T_pre, T_fut = T_fut, S_pre = S_pre, S_fut = S_fut, S_int = S_int,
    D_min = D_min, D_eq = D_eq, D_eq_pre = D_eq_pre,
    H = H, H_med = H_med, W = W, c_lo_banda = cc[1], c_canon = cc[2], c_hi_banda = cc[3],
    c_pkg_lo = c_lo, c_pkg_hi = c_hi,
    A_lo = A_th[1], A_mid = A_th[2], A_hi = A_th[3],
    S_fut_lo = S_th[1], S_fut_mid = S_th[2], S_fut_hi = S_th[3],
    int_estable_km2 = int_es, int_perdido_km2 = int_pe, int_nuevo_km2 = int_nu,
    region = aps$region[i], region_pre_km2 = reg_pre, region_fut_km2 = reg_fut, region_estable_km2 = reg_est,
    chile_pre_km2 = a_pre, chile_fut_km2 = a_fut, ext_km2 = a_ext,
    A_bajo_min = as.integer(n_A < A_MINIMO),
    stringsAsFactors = FALSE)
  list(fila = fila, matriz = M_b, matriz_th = M_th)
}

res <- lapply(seq_len(nrow(aps)), procesar_ap)
completa    <- do.call(rbind, lapply(res, `[[`, "fila"))
matrices    <- do.call(rbind, lapply(res, `[[`, "matriz"))
matrices_th <- do.call(rbind, lapply(res, `[[`, "matriz_th"))
message("Curvas listas | ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")


# 5 Controles ----

message("\n--- Controles ---")
tol <- 1e-5   # los rasteres de Mahalanobis se guardan en FLT4S
dentro <- with(completa, c_canon >= c_pkg_lo * (1 - tol) & c_canon <= c_pkg_hi * (1 + tol))
print(completa[, c("ap", "c_pkg_lo", "c_canon", "c_pkg_hi")], digits = 7)
if (!all(dentro))
  stop("El radio reconstruido no reproduce la clasificación de mh_rep_ch. Revisar la versión de la función o el polígono.")
if (any(abs(completa$A_mid / completa$A - 1) > 1e-3))
  warning("La curva central de la banda no reproduce A. Revisar la coherencia entre Change y Mh_Raw.")
message("Fracción A / celdas interiores (debe rondar ", TH_CH, "): ",
        paste(completa$ap, round(completa$frac_A, 3), collapse = " | "))
if (any(completa$A_bajo_min == 1))
  message("Envolvente estimada con menos de ", A_MINIMO, " celdas en ",
          paste(completa$ap[completa$A_bajo_min == 1], collapse = ", "),
          ". El sesgo va hacia más aislamiento; esas cifras se presentan con reserva.")
vuelca <- with(completa, (S_fut_lo < 1) != (S_fut_hi < 1))
if (any(vuelca, na.rm = TRUE))
  message("La bandera cambia de estado dentro de la banda en ", paste(completa$ap[vuelca], collapse = ", "))


# 6 Tablas y frases ----

indices <- with(completa, data.frame(
  espacio = nombre, ap = ap, area_km2 = round(area_km2, 1),
  A_km2 = ifelse(A_bajo_min == 1, paste0(round(A, 1), " *"), as.character(round(A, 1))),
  T_pre_pct = signif(100 * T_pre, 3), T_fut_pct = signif(100 * T_fut, 3),
  H = signif(H, 2), W = round(W, 1), S_int = round(S_int, 2),
  D_min_km = round(D_min, 1), D_eq_km = round(D_eq, 1), S_fut = round(S_fut, 2),
  dT = round(dT, 3), L = round(L, 1), bandera = ifelse(bandera == 1, "SI", ""),
  escenario = ETIQ_ESC, umbral = TH_CH))

gestion <- with(completa, data.frame(
  espacio = nombre, ap = ap,
  interior_conserva_pct = round(100 * int_estable_km2 / A, 1),
  interior_pierde_pct   = round(100 * int_perdido_km2 / A, 1),
  region = region, region_presente_km2 = round(region_pre_km2), region_futuro_km2 = round(region_fut_km2),
  region_conserva_pct = round(100 * region_estable_km2 / region_pre_km2, 1),
  chile_presente_km2 = round(chile_pre_km2), chile_futuro_km2 = round(chile_fut_km2)))

frases <- with(completa, data.frame(espacio = nombre, ap = ap, frase = ifelse(
  bandera == 1,
  sprintf(paste("El clima actual de %s se da en el presente en el %.2f %% del territorio chileno exterior al área.",
                "En el escenario %s no queda en Chile una superficie con ese clima equivalente al tamaño",
                "del área. Dentro del área conserva su clima el %.0f %% de la superficie."),
          nombre, 100 * T_pre, ETIQ_ESC, 100 * int_estable_km2 / A),
  sprintf(paste("El clima actual de %s se da en el presente en el %.2f %% del territorio chileno exterior al área y",
                "pasa al %.2f %% en el escenario %s. Dentro del área conserva su clima el %.0f %% de la",
                "superficie. Una superficie equivalente al área aparece a %.0f km de su límite, %.1f veces",
                "más lejos de lo que permite la geometría. En la región %s, la superficie con ese clima",
                "pasa de %.0f a %.0f km2."),
          nombre, 100 * T_pre, 100 * T_fut, ETIQ_ESC, 100 * int_estable_km2 / A, D_eq, L,
          region, region_pre_km2, region_fut_km2))))

message("\n=== Índices de aislamiento ===");  print(indices)
message("\n=== Resumen para gestión ===");    print(gestion)
message("\n=== Frases ===");  for (k in seq_len(nrow(frases))) message("- ", frases$frase[k])


# 7 Figuras ----

ESCALA <- c(0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50,
            100, 200, 500, 1000, 2000, 5000, 10000)
fmt <- function(x) format(x, scientific = FALSE, trim = TRUE, drop0trailing = TRUE)
COL_DISP <- "grey45"; COL_PRE <- "#1F6FB4"; COL_FUT <- "#C8501E"
COL_BANDA <- "#E8B9A538"; COL_BANDA2 <- "#C8501E55"
SEP_MIN_REL <- 0.10

ajusta_ancho <- function(txt, cex, ancho_in) {
  pal <- unlist(strsplit(txt, " ", fixed = TRUE)); pal <- pal[nzchar(pal)]
  lin <- character(0); act <- ""
  for (w in pal) {
    pr <- if (nchar(act)) paste(act, w) else w
    if (nchar(act) && strwidth(pr, units = "inches", cex = cex) > ancho_in) { lin <- c(lin, act); act <- w } else act <- pr
  }
  c(lin, act)
}

curva_de <- function(codigo) {
  f <- completa[completa$ap == codigo, ]
  g <- matrices[matrices$ap == codigo, ]; g <- g[order(g$banda), ]
  cur <- data.frame(d = g$d_km, disp = cumsum(g$disponibles) / f$A,
                    pre = cumsum(g$stable + g$lost) / f$A, fut = cumsum(g$stable + g$novel) / f$A)
  cur <- cur[cur$d > 0 & cur$disp > 0, ]
  cur$pre[cur$pre == 0] <- NA; cur$fut[cur$fut == 0] <- NA
  cur
}

# 7.1 Curva de cada AP, autodescriptiva como en el método de aislamiento.
dibuja_curva <- function(codigo) {
  f <- completa[completa$ap == codigo, ]
  cur <- curva_de(codigo)
  if (nrow(cur) < 2) return(invisible(NULL))
  gb <- matrices_th[matrices_th$ap == codigo, ]; gb <- gb[order(gb$banda), ]
  banda <- data.frame(d = gb$d_km, lo = cumsum(gb$fut_lo) / f$A_lo, hi = cumsum(gb$fut_hi) / f$A_hi)
  banda <- banda[banda$d > 0 & banda$lo > 0 & banda$hi > 0, ]
  if (nrow(banda) < 2) banda <- NULL

  ymin <- max(1e-3, min(cur$disp)); ymax <- max(c(cur$disp, 1, banda$hi)) * 1.4
  xmin <- min(cur$d); xmax <- max(cur$d)
  hy <- function(p) ymin * (ymax / ymin)^p
  ancho_log <- log10(xmax / xmin)
  sep_log <- if (is.finite(f$D_min) && is.finite(f$D_eq) && f$D_min > 0) abs(log10(f$D_eq / f$D_min)) else Inf
  juntos <- is.finite(sep_log) && sep_log < SEP_MIN_REL * ancho_log

  png(file.path(DIR_FIG, paste0("curva_", codigo, ".png")), width = 1700, height = 1200, res = 150)
  on.exit(dev.off(), add = TRUE)
  par(mar = c(8.8, 4.8, 5.8, 1.4), mgp = c(2.8, 0.7, 0))
  plot(NA, log = "xy", xlim = c(xmin, xmax), ylim = c(ymin, ymax), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  ty <- ESCALA[ESCALA >= ymin & ESCALA <= ymax]; tx <- ESCALA[ESCALA >= xmin & ESCALA <= xmax]
  axis(2, at = ty, labels = fmt(ty), las = 1, cex.axis = 0.8, tcl = -0.3)
  axis(2, at = 1, labels = "1", las = 1, font = 2, cex.axis = 0.9, tcl = -0.45)
  axis(1, at = tx, labels = fmt(tx), cex.axis = 0.8, tcl = -0.3)
  mtext("Distancia al límite del área protegida (km)", side = 1, line = 2.5, cex = 0.95)
  mtext("Superficie acumulada (múltiplos de A)", side = 2, line = 3.2, cex = 0.92)

  mtext(sprintf("%s   [%s]", f$nombre, codigo), side = 3, line = 4.3, font = 2, cex = 1.08)
  mtext(paste("Clima futuro", ETIQ_ESC), side = 3, line = 3.1, cex = 0.8, col = "grey30")
  sub1 <- if (f$bandera == 1)
    sprintf("dT = %.2f    SIN SUSTITUTO CLIMÁTICO del tamaño del área (S_fut = %.2f < 1)", f$dT, f$S_fut) else
      sprintf("dT = %.2f    L = %.1f    S_fut = %.2f", f$dT, f$L, f$S_fut)
  mtext(sub1, side = 3, line = 1.8, cex = 0.88, col = if (f$bandera == 1) "#B03A12" else "grey20")
  mtext(sprintf("A = %.1f km2 (%d celdas%s)    T_pre = %s %% del territorio exterior    S_int = %.2f    D_min = %.1f km    D_eq = %s    H = %s    W = %.1f",
                f$A, f$n_A, if (f$A_bajo_min == 1) ", envolvente con pocos datos" else "",
                fmt(signif(100 * f$T_pre, 3)), f$S_int, f$D_min,
                if (is.finite(f$D_eq)) sprintf("%.1f km", f$D_eq) else "no se alcanza",
                fmt(signif(f$H, 2)), f$W),
        side = 3, line = 0.5, cex = 0.7, col = "grey35")

  if (!is.null(banda)) {
    polygon(c(banda$d, rev(banda$d)), c(pmax(banda$hi, ymin), rev(pmax(banda$lo, ymin))), col = COL_BANDA, border = NA)
    lines(banda$d, pmax(banda$lo, ymin), col = COL_BANDA2, lwd = 0.9, lty = 3)
    lines(banda$d, pmax(banda$hi, ymin), col = COL_BANDA2, lwd = 0.9, lty = 3)
  }
  abline(h = 1, col = "grey55", lty = 2)
  text(xmax, 1, "A", adj = c(1.2, -0.5), cex = 0.85, font = 2, col = "grey35")

  # dT entre mesetas, anclado donde ambas curvas ya son meseta.
  if (f$S_pre > 0 && f$S_fut > 0) {
    ini <- function(y) { y <- ifelse(is.na(y), 0, y); fin <- y[length(y)]
    if (fin <= 0) NA_real_ else cur$d[min(which(y >= 0.99 * fin))] }
    d_ini <- suppressWarnings(max(ini(cur$pre), ini(cur$fut)))
    if (!is.finite(d_ini) || d_ini <= xmin || log10(xmax / d_ini) < 0.05 * ancho_log)
      d_ini <- 10^(log10(xmax) - 0.10 * ancho_log)
    y1 <- min(f$S_pre, f$S_fut); y2 <- max(f$S_pre, f$S_fut)
    rect(d_ini, y1, xmax, y2, col = "#1F6FB433", border = NA)
    text(10^mean(log10(c(d_ini, xmax))), 10^mean(log10(c(y1, y2))), "dT", cex = 0.95, font = 2, col = COL_PRE)
  }
  if (is.finite(f$D_min) && is.finite(f$D_eq) && f$D_eq > f$D_min) {
    rect(f$D_min, ymin, f$D_eq, 1, col = "#F2C3B055", border = NA)
    if (!juntos) text(10^mean(log10(c(f$D_min, f$D_eq))), hy(0.15), "L", cex = 0.95, font = 2, col = "#8C3A14")
  }
  lines(cur$d, cur$disp, col = COL_DISP, lwd = 2, lty = 3)
  lines(cur$d, cur$pre, col = COL_PRE, lwd = 2.2)
  lines(cur$d, cur$fut, col = COL_FUT, lwd = 2.6)
  for (s in list(list(v = f$S_pre, c = COL_PRE, t = "S_pre"), list(v = f$S_fut, c = COL_FUT, t = "S_fut")))
    if (s$v >= ymin && s$v <= ymax) {
      abline(h = s$v, col = s$c, lty = 3)
      text(xmin, s$v, sprintf("%s = %.2f", s$t, s$v), adj = c(-0.05, -0.45), cex = 0.72, col = s$c)
    }
  if (is.finite(f$D_min)) {
    segments(f$D_min, ymin, f$D_min, 1, col = COL_DISP, lty = 3)
    points(f$D_min, 1, pch = 21, bg = COL_DISP, col = "white", cex = 1.3)
    text(f$D_min, hy(0.04), expression(D[min]), adj = if (juntos) c(1.15, 0.5) else c(-0.15, 0.5),
         cex = 0.85, font = 2, col = "grey30")
  }
  if (is.finite(f$D_eq)) {
    segments(f$D_eq, ymin, f$D_eq, 1, col = COL_FUT, lty = 2)
    points(f$D_eq, 1, pch = 21, bg = COL_FUT, col = "white", cex = 1.3)
    text(f$D_eq, if (juntos) hy(0.13) else hy(0.04), expression(D[eq]), adj = c(-0.15, 0.5),
         cex = 0.85, font = 2, col = COL_FUT)
  } else if (f$bandera == 1) {
    text(xmax, hy(0.55), "la curva de análogo futuro no alcanza A", adj = c(1, 0.5), cex = 0.82, font = 2, col = "#B03A12")
  }
  if (is.finite(f$D_eq_pre)) points(f$D_eq_pre, 1, pch = 21, bg = "white", col = COL_PRE, cex = 1.1, lwd = 1.6)
  d0 <- suppressWarnings(min(cur$d[is.finite(cur$fut)]))
  if (is.finite(d0) && d0 > xmin * 1.5)
    text(xmin, hy(0.30), sprintf("sin análogo futuro a menos de %.0f km del límite", d0),
         adj = c(-0.02, 0.5), cex = 0.74, font = 3, col = COL_FUT)

  leg <- c("Territorio disponible (cualquier clima)", "Análogo presente", "Análogo futuro")
  cols <- c(COL_DISP, COL_PRE, COL_FUT); ltys <- c(3, 1, 1); lwds <- c(2, 2.2, 2.6)
  if (!is.null(banda)) {
    leg <- c(leg, sprintf("Análogo futuro con envolvente a %s y %s (suplementario)", TH_BANDA[1], TH_BANDA[3]))
    cols <- c(cols, COL_BANDA2); ltys <- c(ltys, 3); lwds <- c(lwds, 1.2)
  }
  legend("topleft", bty = "o", bg = "white", box.col = "grey85", cex = 0.76, lwd = lwds, lty = ltys,
         seg.len = 2.4, col = cols, legend = leg)

  pie <- c(paste("Cada curva acumula la superficie al alejarse del límite del área, en múltiplos de A (superficie",
                 "interior con su propio clima). Los puntos llenos marcan dónde cada curva alcanza A y el hueco,",
                 "el análogo presente. Sombreado horizontal L (D_eq / D_min), vertical dT (S_fut / S_pre)."),
           if (!is.null(banda)) sprintf(paste("Banda: la misma curva futura con el radio del envolvente en los cuantiles",
                                              "%s y %s de la D2 interior (c = %s, %s, %s). Cada curva va escalada por su propio A.",
                                              "Es sensibilidad a la definición, no un intervalo de confianza."),
                                        TH_BANDA[1], TH_BANDA[3], fmt(signif(f$c_lo_banda, 3)),
                                        fmt(signif(f$c_canon, 3)), fmt(signif(f$c_hi_banda, 3))),
           sprintf(paste("Clima futuro %s. Índices de la clasificación canónica al cuantil %s.",
                         "Superficies ponderadas por el área real de cada celda."), ETIQ_ESC, TH_CH))
  lp <- unlist(lapply(pie, ajusta_ancho, cex = 0.64, ancho_in = par("pin")[1]))[1:7]
  lp <- lp[!is.na(lp)]
  for (j in seq_along(lp)) mtext(lp[j], side = 1, line = 3.9 + 0.8 * (j - 1), cex = 0.64, col = "grey40", adj = 0)
  box(col = "grey60")
  invisible(NULL)
}
for (k in completa$ap) dibuja_curva(k)

# 7.2 Comparación de las AP en una figura, para la diapositiva. Presente
# discontinuo y futuro continuo, cada AP escalada por su propio A.
COL_AP <- c("cadetblue", "darkorchid", "brown1", "chartreuse3")
png(file.path(DIR_FIG, "comparacion_curvas_Aysen.png"), width = 1700, height = 1150, res = 150)
par(mar = c(5, 5, 4.5, 1.5))
curvas <- lapply(completa$ap, curva_de)
xr <- range(unlist(lapply(curvas, `[[`, "d")))
yr <- range(c(unlist(lapply(curvas, function(x) c(x$pre, x$fut))), 1), na.rm = TRUE) * c(0.8, 1.4)
plot(NA, log = "xy", xlim = xr, ylim = yr, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
ty <- ESCALA[ESCALA >= yr[1] & ESCALA <= yr[2]]; tx <- ESCALA[ESCALA >= xr[1] & ESCALA <= xr[2]]
axis(2, at = ty, labels = fmt(ty), las = 1, cex.axis = 0.85); axis(1, at = tx, labels = fmt(tx), cex.axis = 0.85)
mtext("Distancia al límite del área protegida (km)", side = 1, line = 2.8)
mtext("Superficie con clima equivalente (veces el tamaño del área)", side = 2, line = 3.4)
mtext("Superficie con el clima actual de cada área, presente y en el futuro", side = 3, line = 2.4, font = 2, cex = 1.1)
mtext(paste("Clima futuro", ETIQ_ESC), side = 3, line = 1.1, cex = 0.85, col = "grey30")
abline(h = 1, col = "grey55", lty = 2)
for (i in seq_along(curvas)) {
  lines(curvas[[i]]$d, curvas[[i]]$pre, col = COL_AP[i], lwd = 2, lty = 2)
  lines(curvas[[i]]$d, curvas[[i]]$fut, col = COL_AP[i], lwd = 2.8)
  if (is.finite(completa$D_eq[i])) points(completa$D_eq[i], 1, pch = 21, bg = COL_AP[i], col = "white", cex = 1.5)
}
legend("topleft", bty = "o", bg = "white", box.col = "grey85", cex = 0.85,
       legend = c(sprintf("%s, presente", completa$nombre), sprintf("%s, futuro", completa$nombre)),
       col = rep(COL_AP[seq_along(curvas)], 2), lty = rep(c(2, 1), each = length(curvas)), lwd = rep(c(2, 2.8), each = length(curvas)))
dev.off()

# 7.3 Mapas de cambio, uno por AP y centrado en ella. La ventana es un cuadrado
# en LAEA que contiene la AP con un margen de MARGEN_MAPA_KM, de modo que todos
# los paneles tienen la misma forma aunque las AP difieran en tamaño. El mapa de
# todo Chile de cada AP lo genera mh_rep_ch en Charts.
COLS_CAMBIO <- c("grey90", "aquamarine4", "coral1", "steelblue2")
ETIQ_CAMBIO <- c("No análogo", "Estable", "Perdido", "Nuevo")
ventana_ap <- function(i) {
  bb <- st_bbox(st_buffer(st_transform(st_geometry(aps)[i], CRS_EA), MARGEN_MAPA_KM * 1000))
  cx <- (bb[["xmin"]] + bb[["xmax"]]) / 2; cy <- (bb[["ymin"]] + bb[["ymax"]]) / 2
  h  <- max(bb[["xmax"]] - bb[["xmin"]], bb[["ymax"]] - bb[["ymin"]]) / 2
  sq <- st_as_sfc(st_bbox(c(xmin = cx - h, ymin = cy - h, xmax = cx + h, ymax = cy + h), crs = st_crs(CRS_EA)))
  ext(vect(st_sf(geometry = st_transform(sq, 4326))))
}
n_col <- ceiling(sqrt(nrow(aps))); n_fil <- ceiling(nrow(aps) / n_col)
png(file.path(DIR_FIG, "mapas_cambio_Aysen.png"), width = 900 * n_col, height = 900 * n_fil + 180, res = 150)
par(mfrow = c(n_fil, n_col), oma = c(3, 0, 2.5, 0))
for (i in seq_len(nrow(aps))) {
  ven <- ventana_ap(i)
  rz  <- crop(rast(f_ch[i]), ven)
  plot(rz, breaks = c(-0.5, 0.5, 1.5, 2.5, 3.5), col = COLS_CAMBIO, legend = FALSE,
       main = sprintf("%s [%s]", aps$nombre[i], aps[[COL_ID]][i]))
  lines(vect(reg_ap[[i]]), col = "grey25", lwd = 0.8)
  if (nrow(aps) > 1) lines(vect(aps[-i, ]), col = "grey35", lty = 2)
  lines(vect(aps[i, ]), col = "black", lwd = 2)
}
mtext(paste("Clima análogo al de cada área. Clima futuro", ETIQ_ESC), outer = TRUE, side = 3, line = 0.8, font = 2)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend("bottom", horiz = TRUE, fill = COLS_CAMBIO, legend = ETIQ_CAMBIO, bty = "n", cex = 0.95)
dev.off()


# 8 Salidas ----

write.csv(indices,     file.path(DIR_CH, "indices_aislamiento_Aysen.csv"), row.names = FALSE)
write.csv(gestion,     file.path(DIR_CH, "resumen_gestion_Aysen.csv"),     row.names = FALSE)
write.csv(frases,      file.path(DIR_CH, "frases_gestion_Aysen.csv"),      row.names = FALSE)
write.csv(completa,    file.path(DIR_CH, "tabla_completa_Aysen.csv"),      row.names = FALSE)
write.csv(matrices,    file.path(DIR_CH, "matrices_bandas_Aysen.csv"),     row.names = FALSE)
write.csv(matrices_th, file.path(DIR_CH, "matrices_bandas_th_Aysen.csv"),  row.names = FALSE)
write_sf(merge(aps[, COL_ID], completa, by.x = COL_ID, by.y = "ap"), file.path(DIR_CH, "APs_Aysen_aislamiento.gpkg"),
         delete_dsn = TRUE)
writeLines(capture.output(sessionInfo()), file.path(DIR_CH, "sessionInfo.txt"))
message("\nListo en ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min | ", DIR_CH)
