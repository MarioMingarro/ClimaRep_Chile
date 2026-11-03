suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(ClimaRep)
  library(parallel)
  library(foreach)
  library(doParallel)
})

# Representatividad climática de las AP de Chile por superficie.
# Fondo (BCRS) y AP (APCRS) se construyen con el mismo estimador. Malla de
# UNIT_KM en LAEA; cada cuadrícula del fondo es la tierra que contiene y cada
# fragmento protegido es la red de AP intersecada con la cuadrícula y disuelta.
# Cada unidad o fragmento tiene su envolvente P99 (regla de celdas de ClimaRep)
# y suma su superficie / UNIT_KM^2 en las celdas de su envolvente.
# RCRI por celda en la rejilla equivalente de 1 km; huérfanos fuera del RCRI y
# en una capa propia.


# 0 Configuración ----

DIR_BASE        <- "C:/A_TRABAJO/CHILE"
DIR_CHELSA      <- "C:/A_TRABAJO/DATA/CHELSA/BIOCLIM_1981_2010"
PATH_STUDY_AREA <- file.path(DIR_BASE, "DATA/Limites_admin/chile_limites.shp")
PATH_AP_RAW     <- file.path(DIR_BASE, "DATA/Areas_Protegidas/APs_chile_2026.shp")
PATH_AP         <- file.path(DIR_BASE, "DATA/Areas_Protegidas/APs_chile_2026_5km2.shp")

COL_ID          <- "SITECODE"
TH              <- 0.99       # umbral de las envolventes (mh_rep y BCRS)
VIF_TH          <- 5
MIN_AREA_M2     <- 5e6        # 5 km2, tras recortar con el área de estudio
UNIT_KM         <- 10         # lado de la malla de referencia
RES_EA          <- 1000       # resolución de la rejilla equivalente final (m)
NUM_CORES       <- 8
RCRI_LOG        <- log        # cambiar a log10 si así se definió en Natura 2000

RUN_REPLICATION <- TRUE       # versión complementaria por AP (cada AP suma 1)
N_VAL_RANDOM    <- 30         # validación contra mh_rep(), fragmentos al azar
N_VAL_SMALL     <- 15         # fragmentos más pequeños con celdas
N_VAL_MULTI     <- 10         # fragmentos multiparte

# LAEA centrada en Chile continental, elipsoide GRS80 (SIRGAS).
CRS_EA <- "+proj=laea +lat_0=-37 +lon_0=-71 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"

P        <- round(TH * 100)
DIR_OUT  <- file.path(DIR_BASE, "RESULTS", paste0("APCRS_P", P))
DIR_RUN  <- file.path(DIR_OUT, "mh_rep")
DIR_TMP  <- file.path(DIR_OUT, "tmp")
CLIM_VIF <- file.path(DIR_OUT, "present_climate_VIF_WGS84.tif")
SA_GPKG  <- file.path(DIR_OUT, "study_area_WGS84.gpkg")
for (d in c(DIR_RUN, DIR_TMP)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
terraOptions(tempdir = DIR_TMP)
sf_use_s2(FALSE)
t0  <- Sys.time()
s_m <- UNIT_KM * 1000

TAG_BA <- sprintf("BCRS_area_%dkm_P%d", UNIT_KM, P)    # fondo por superficie (principal)
TAG_AA <- sprintf("APCRS_area_%dkm_P%d", UNIT_KM, P)   # AP por superficie (principal)
TAG_BI <- sprintf("BCRS_%dkm_P%d", UNIT_KM, P)         # fondo entero (complemento)


# 1 Área de estudio y rejilla equivalente ----
# La máscara exacta se aplica una vez al clima. A mh_rep() se le pasa el
# rectángulo envolvente, que da el mismo espacio de referencia porque el
# clima ya es NA fuera de la máscara.

sa <- read_sf(PATH_STUDY_AREA) |> st_make_valid() |> st_transform(4326) |> st_union() |> st_make_valid()
sa <- st_sf(ID = 1L, geometry = sa)
sa_box <- st_sf(ID = 1L, geometry = st_as_sfc(st_bbox(sa)))
write_sf(sa_box, SA_GPKG, delete_dsn = file.exists(SA_GPKG))
sa_ea_sf <- st_make_valid(st_transform(sa, CRS_EA))
sa_ea    <- vect(sa_ea_sf)

e <- as.vector(ext(sa_ea))
tpl <- rast(xmin = floor(e[1] / RES_EA) * RES_EA, xmax = ceiling(e[2] / RES_EA) * RES_EA,
            ymin = floor(e[3] / RES_EA) * RES_EA, ymax = ceiling(e[4] / RES_EA) * RES_EA,
            resolution = RES_EA, crs = CRS_EA)
writeRaster(init(tpl, 0L), file.path(DIR_OUT, "template_EA_1km.tif"), datatype = "INT1U",
            overwrite = TRUE, gdal = "COMPRESS=LZW")


# 2 Clima ----
# BIO1 temperatura media anual; BIO2 rango diurno medio; BIO4 estacionalidad
# de la temperatura (desviación típica x 100); BIO12 precipitación anual;
# BIO15 estacionalidad de la precipitación (coeficiente de variación).
# Variables sin transformar.

VARS <- c("bio1", "bio2", "bio4", "bio12", "bio15")

if (!file.exists(CLIM_VIF)) {
  f <- list.files(DIR_CHELSA, "bio[0-9]+.*\\.tif$", full.names = TRUE, ignore.case = TRUE)
  v <- sub("^bio0+([1-9])", "bio\\1", sub("^.*?(bio[0-9]+).*$", "\\1", tolower(basename(f)), perl = TRUE))
  clim <- rast(f[match(VARS, v)]); names(clim) <- VARS
  crs(clim) <- "EPSG:4326"
  sa_v <- vect(sa)
  clim <- mask(crop(clim, sa_v, snap = "out"), sa_v, touches = FALSE)
  vif <- vif_filter(clim, th = VIF_TH)                 # solo diagnóstico
  message("VIF retendría: ", paste(vif$summary$kept_layers, collapse = ", "))
  writeRaster(clim, CLIM_VIF, overwrite = TRUE,
              gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "BIGTIFF=IF_SAFER"))
}
message("Variables: ", paste(names(rast(CLIM_VIF)), collapse = ", "))


# 3 Áreas protegidas ----
# Una geometría por COL_ID, recorte con el área de estudio y filtro de área
# en la proyección equivalente. Sin filtro de IPQ.

aps <- read_sf(PATH_AP_RAW) |> st_make_valid() |> st_transform(CRS_EA) |> st_make_valid()
aps <- aggregate(aps[, COL_ID], by = list(id = aps[[COL_ID]]), FUN = function(x) x[1])[, COL_ID]
aps <- st_intersection(aps, st_geometry(sa_ea_sf))
aps <- st_collection_extract(aps, "POLYGON", warn = FALSE)
aps <- aggregate(aps[, COL_ID], by = list(id = aps[[COL_ID]]), FUN = function(x) x[1])[, COL_ID]
aps <- aps[as.numeric(st_area(aps)) > MIN_AREA_M2, ]
aps_ea <- aps
aps <- st_transform(aps, 4326)                       # mh_rep() exige EPSG:4326
write_sf(aps, PATH_AP, delete_dsn = file.exists(PATH_AP))
if (anyDuplicated(tolower(aps[[COL_ID]]))) stop("Hay ", COL_ID, " duplicados.")
message("Áreas protegidas > 5 km2: ", nrow(aps))


# 4 Clima blanqueado y asignación de celdas ----
# Cada celda pertenece a la cuadrícula que contiene su centro y está protegida
# si su centro cae dentro de alguna AP (regla de mask(crop()) de ClimaRep).
# D2 se calcula en espacio blanqueado con la covarianza de toda el área de
# estudio y se poda por el primer eje (poda exacta).

clim <- rast(CLIM_VIF)
X  <- values(clim, mat = TRUE)
ok <- which(complete.cases(X))
X  <- X[ok, , drop = FALSE]
Z  <- X %*% solve(chol(cov(X)))
rm(X); invisible(gc())

xy_wgs <- xyFromCell(clim, ok)
xy_ea  <- project(xy_wgs, from = crs(clim), to = CRS_EA)
uid    <- paste(floor(xy_ea[, 1] / s_m), floor(xy_ea[, 2] / s_m), sep = "_")

prot_r <- rasterize(vect(aps), clim[[1]], field = 1, background = 0, touches = FALSE)
prot   <- values(prot_r, mat = FALSE)[ok] == 1
rm(prot_r)

o   <- order(Z[, 1])
Zs  <- Z[o, , drop = FALSE]; rm(Z)
z1s <- Zs[, 1]
z2s <- rowSums(Zs^2)
pos <- integer(length(o)); pos[o] <- seq_along(o)
idx_u <- split(pos, uid)                  # cuadrículas del fondo
idx_f <- split(pos[prot], uid[prot])      # fragmentos protegidos
message(sprintf("Celdas con clima %d | protegidas %d | cuadrículas %d | con celdas protegidas %d",
                length(ok), sum(prot), length(idx_u), length(idx_f)))


# 5 Malla y pesos geométricos ----
# Pesos a ambos lados con la misma definición de superficie, el área en LAEA
# dividida entre UNIT_KM^2. Solo se consideran cuadrículas con alguna celda de
# clima. Las cuadrículas contenidas por completo valen 1 y solo se intersecan
# las de borde.

sq_from_uid <- function(u) {
  ij <- matrix(as.numeric(unlist(strsplit(u, "_"))), ncol = 2, byrow = TRUE)
  st_sf(uid = u, geometry = st_sfc(lapply(seq_len(nrow(ij)), function(i) {
    x0 <- ij[i, 1] * s_m; y0 <- ij[i, 2] * s_m
    st_polygon(list(rbind(c(x0, y0), c(x0 + s_m, y0), c(x0 + s_m, y0 + s_m),
                          c(x0, y0 + s_m), c(x0, y0))))
  }), crs = CRS_EA))
}

units_u <- names(idx_u)
sq_u    <- sq_from_uid(units_u)

# Tierra de cada cuadrícula
g_sa  <- st_geometry(sa_ea_sf)
hit   <- lengths(st_intersects(sq_u, g_sa)) > 0
full  <- hit & lengths(st_within(sq_u, g_sa)) > 0
land  <- setNames(numeric(nrow(sq_u)), units_u)
land[full] <- s_m^2
brd <- which(hit & !full)
if (length(brd)) {
  xl <- suppressWarnings(st_intersection(sq_u[brd, ], g_sa))
  al <- tapply(as.numeric(st_area(xl)), xl$uid, sum)
  land[names(al)] <- al
}
w_u <- pmin(land / s_m^2, 1)
if (any(w_u == 0)) message("Cuadrículas con celdas de clima y sin tierra en la geometría: ", sum(w_u == 0))

# Fragmentos protegidos, intersección de las AP con la malla y disolución por
# cuadrícula (elimina solapes entre AP)
xp <- suppressWarnings(st_intersection(sq_u[, "uid"], st_geometry(aps_ea)))
xp <- st_collection_extract(st_make_valid(xp), "POLYGON", warn = FALSE)
frag <- aggregate(xp[, "uid"], by = list(g = xp$uid), FUN = function(x) x[1])[, "uid"]
stopifnot(!anyDuplicated(frag$uid))
frag_a <- setNames(as.numeric(st_area(frag)), frag$uid)
w_f    <- pmin(frag_a / s_m^2, 1)
parts  <- suppressWarnings(st_cast(st_cast(frag, "MULTIPOLYGON"), "POLYGON"))
n_parts <- table(parts$uid)
rm(xp)

# Coherencia entre geometría y celdas
sin_geom  <- setdiff(names(idx_f), names(w_f))
sin_celda <- setdiff(names(w_f), names(idx_f))
if (length(sin_geom)) message("Cuadrículas con celdas protegidas y sin fragmento geométrico: ",
                              length(sin_geom), " (", sum(lengths(idx_f[sin_geom])), " celdas, se descartan)")
frags <- intersect(names(idx_f), names(w_f))
a_ap_total <- sum(as.numeric(st_area(aps_ea)))
message(sprintf(paste("Fragmentos %d | multiparte %d | sin ningún centro de celda %d (%.2f km2)",
                      "| superficie de fragmentos %.0f km2 | suma de áreas de AP %.0f km2",
                      "| fracción protegida %.4f"),
                length(frags), sum(n_parts[frags] > 1), length(sin_celda),
                sum(frag_a[sin_celda]) / 1e6, sum(frag_a) / 1e6, a_ap_total / 1e6,
                sum(w_f[frags]) / sum(w_u)))


# 6 Envolventes, BCRS y APCRS ----

env_hits <- function(p) {
  zu  <- Zs[p, , drop = FALSE]
  zm  <- colMeans(zu)
  d2u <- rowSums(sweep(zu, 2, zm)^2)
  thr <- if (TH >= 1) max(d2u) else as.numeric(quantile(d2u, TH, type = 7))
  thr <- thr * (1 + 1e-9) + 1e-12
  lo  <- findInterval(zm[1] - sqrt(thr), z1s) + 1L
  hi  <- findInterval(zm[1] + sqrt(thr), z1s)
  if (hi < lo) return(list(hit = integer(0), thr = thr))
  j  <- lo:hi
  d2 <- z2s[j] - 2 * drop(Zs[j, , drop = FALSE] %*% zm) + sum(zm^2)
  list(hit = j[d2 <= thr], thr = thr)
}

bcrs_a  <- numeric(length(z1s))   # fondo por superficie
bcrs_i  <- integer(length(z1s))   # fondo entero (complemento de replicación)
apcrs_a <- numeric(length(z1s))   # AP por superficie

for (k in seq_along(units_u)) {
  h <- env_hits(idx_u[[units_u[k]]])$hit
  bcrs_i[h] <- bcrs_i[h] + 1L
  bcrs_a[h] <- bcrs_a[h] + w_u[[units_u[k]]]
  if (k %% 1000 == 0) message("BCRS ", k, " / ", length(units_u))
}

thr_f <- setNames(numeric(length(frags)), frags)
for (k in seq_along(frags)) {
  e_k <- env_hits(idx_f[[frags[k]]])
  apcrs_a[e_k$hit] <- apcrs_a[e_k$hit] + w_f[[frags[k]]]
  thr_f[k] <- e_k$thr
  if (k %% 500 == 0) message("APCRS ", k, " / ", length(frags))
}

write_vec <- function(vs, tag, dtype) {
  r <- rast(clim, nlyrs = 1); names(r) <- "value"
  vals <- rep(NA, ncell(r)); vals[ok] <- vs[pos]; values(r) <- vals
  gd <- if (dtype == "INT4S") "COMPRESS=LZW" else c("COMPRESS=DEFLATE", "PREDICTOR=3")
  f <- file.path(DIR_OUT, paste0(tag, "_WGS84.tif"))
  writeRaster(r, f, datatype = dtype, overwrite = TRUE, gdal = gd)
  r_ea <- mask(project(rast(f), tpl, method = "near"), sa_ea)
  writeRaster(r_ea, file.path(DIR_OUT, paste0(tag, "_EA_1km.tif")),
              datatype = dtype, overwrite = TRUE, gdal = gd)
}
write_vec(bcrs_a,  TAG_BA, "FLT4S")
write_vec(apcrs_a, TAG_AA, "FLT4S")
write_vec(bcrs_i,  TAG_BI, "INT4S")
message(sprintf("Celdas con BCRS = 0 %d | con APCRS por superficie = 0 %d (%.2f %%)",
                sum(bcrs_a == 0), sum(apcrs_a == 0), 100 * mean(apcrs_a == 0)))


# 7 Diagnóstico de fragmentos multiparte ----
# Para cada fragmento con varias partes, distancia máxima entre los centroides
# climáticos de sus partes dividida entre el radio de su envolvente. Valores
# próximos a 2 indican partes con climas en extremos opuestos de la envolvente,
# que entonces abarca climas intermedios no protegidos.

diag <- data.frame(uid = frags, area_km2 = frag_a[frags] / 1e6,
                   n_parts = as.integer(n_parts[frags]), n_celdas = lengths(idx_f[frags]),
                   sep_rel = NA_real_)
multi <- frags[n_parts[frags] > 1]
if (length(multi)) {
  cm  <- which(prot & uid %in% multi)
  pts <- st_as_sf(data.frame(i = cm, x = xy_ea[cm, 1], y = xy_ea[cm, 2]),
                  coords = c("x", "y"), crs = CRS_EA)
  pm  <- parts[parts$uid %in% multi, "uid"]
  pm$part <- seq_len(nrow(pm))
  jn  <- st_drop_geometry(st_join(pts, pm, join = st_intersects, left = FALSE))
  jn  <- jn[jn$uid == uid[jn$i] & !duplicated(jn$i), ]
  zj  <- Zs[pos[jn$i], , drop = FALSE]
  for (u in unique(jn$uid)) {
    s   <- jn$uid == u
    cen <- rowsum(zj[s, , drop = FALSE], jn$part[s]) / as.vector(table(jn$part[s]))
    if (nrow(cen) > 1) diag$sep_rel[diag$uid == u] <- max(dist(cen)) / sqrt(thr_f[[u]])
  }
}
write.csv(diag, file.path(DIR_OUT, sprintf("diagnostico_fragmentos_P%d.csv", P)), row.names = FALSE)
dm <- diag[diag$n_parts > 1, ]
message(sprintf("Multiparte %d de %d | %.2f %% de la superficie protegida | con partes evaluables %d",
                nrow(dm), nrow(diag), 100 * sum(dm$area_km2) / sum(diag$area_km2), sum(!is.na(dm$sep_rel))))
if (any(!is.na(dm$sep_rel))) {
  print(round(quantile(dm$sep_rel, c(0.5, 0.9, 0.99, 1), na.rm = TRUE), 3))
  message(sprintf("Multiparte con sep_rel > 1 %d (%.2f %% de la superficie protegida)",
                  sum(dm$sep_rel > 1, na.rm = TRUE),
                  100 * sum(dm$area_km2[which(dm$sep_rel > 1)]) / sum(diag$area_km2)))
}


# 8 Validación contra mh_rep() en una muestra de fragmentos ----
# Al azar, los más pequeños y multiparte. Los fragmentos se densifican en LAEA
# antes de pasar a 4326 para que los lados de la malla sigan siendo rectos.

set.seed(1)
small <- frags[order(frag_a[frags])][seq_len(min(N_VAL_SMALL, length(frags)))]
mult  <- if (length(multi)) sample(multi, min(N_VAL_MULTI, length(multi))) else character(0)
rand  <- sample(setdiff(frags, c(small, mult)), min(N_VAL_RANDOM, length(frags)))
s_val <- unique(c(rand, small, mult))

fv <- frag[frag$uid %in% s_val, ]
fv$fname <- sprintf("f%05d", seq_len(nrow(fv)))
fv4 <- st_transform(st_segmentize(fv, 100), 4326)
d_val <- file.path(DIR_OUT, "validacion_fragmentos")
suppressMessages(ClimaRep::mh_rep(polygon = fv4, col_name = "fname",
                                  climate_variables = rast(CLIM_VIF), study_area = sa_box,
                                  th = TH, dir_output = d_val, save_raw = FALSE))

cmp <- do.call(rbind, lapply(seq_len(nrow(fv)), function(i) {
  u  <- fv$uid[i]
  fr <- file.path(d_val, "Analogues", paste0("Th_rep_", fv$fname[i], ".tif"))
  b  <- integer(length(ok)); b[o[env_hits(idx_f[[u]])$hit]] <- 1L
  if (!file.exists(fr)) {
    a <- rep(NA_real_, length(ok))
  } else {
    r <- rast(fr); a <- r[cellFromXY(r, xy_wgs)][, 1]
  }
  data.frame(uid = u, area_km2 = frag_a[[u]] / 1e6, n_parts = as.integer(n_parts[[u]]),
             n_celdas = length(idx_f[[u]]), analogos_fast = sum(b),
             analogos_mh = sum(a, na.rm = TRUE),
             solo_fast = sum(b == 1 & a == 0, na.rm = TRUE),
             solo_mh = sum(b == 0 & a == 1, na.rm = TRUE), na_mh = sum(is.na(a)))
}))
print(cmp)
write.csv(cmp, file.path(DIR_OUT, sprintf("validacion_fragmentos_mh_rep_P%d.csv", P)), row.names = FALSE)
message(sprintf("Validación | fragmentos %d | con alguna diferencia %d",
                nrow(cmp), sum(cmp$solo_fast > 0 | cmp$solo_mh > 0 | cmp$na_mh > 0)))


# 9 RCRI por superficie ----
# RCRI por celda = log((APCRS / suma APCRS) / (BCRS / suma BCRS)) en la
# rejilla equivalente de 1 km, con sumas simples de celdas y sin constante.
# Estado de cada celda, 0 con RCRI válido, 1 huérfano (APCRS = 0, clima sin
# análogos en la red de AP), 2 con BCRS = 0 y APCRS > 0. El RCRI solo se
# define en el estado 0; los estados 1 y 2 quedan en la capa de estado.

rcri_ea <- function(f_ap, f_b) {
  ap <- rast(f_ap); b <- rast(f_b)
  stopifnot(compareGeom(ap, b), global(is.na(ap) != is.na(b), "sum")[[1]] == 0)
  st <- ifel(ap == 0, 1L, ifel(b == 0, 2L, 0L)); names(st) <- "status"
  r  <- RCRI_LOG((ap / global(ap, "sum", na.rm = TRUE)[[1]]) /
                   (b  / global(b,  "sum", na.rm = TRUE)[[1]]))
  r  <- mask(r, st, maskvalues = c(1, 2)); names(r) <- "rcri"
  list(rcri = r, status = st)
}

write_rcri <- function(res, tag) {
  writeRaster(res$rcri, file.path(DIR_OUT, paste0(tag, "_EA_1km.tif")), datatype = "FLT4S",
              overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
  writeRaster(res$status, file.path(DIR_OUT, paste0(tag, "_status_EA_1km.tif")), datatype = "INT1U",
              overwrite = TRUE, gdal = "COMPRESS=LZW")
  n  <- global(!is.na(res$status), "sum")[[1]]
  nv <- global(res$status == 0, "sum", na.rm = TRUE)[[1]]
  message(sprintf(paste("%s | celdas %d | válidas %d | RCRI > 0 %.1f %% de las válidas",
                        "| huérfanos %d (%.2f %%) | BCRS = 0 %d"),
                  tag, n, nv,
                  100 * global(res$rcri > 0, "sum", na.rm = TRUE)[[1]] / nv,
                  global(res$status == 1, "sum", na.rm = TRUE)[[1]],
                  100 * global(res$status == 1, "sum", na.rm = TRUE)[[1]] / n,
                  global(res$status == 2, "sum", na.rm = TRUE)[[1]]))
}

# Resumen por AP. Celdas por su centro, igual que mask(crop()), sin fracciones.
resumen_ap <- function(res, ap_sf, tag) {
  ap_ea <- st_transform(ap_sf, CRS_EA)
  ex <- terra::extract(c(res$rcri, res$status), vect(ap_ea), ID = TRUE)
  ex$ID <- factor(ex$ID, levels = seq_len(nrow(ap_ea)))
  q <- function(x, p) if (all(is.na(x))) NA_real_ else unname(quantile(x, p, na.rm = TRUE, type = 7))
  out <- do.call(rbind, lapply(split(ex, ex$ID), function(d) data.frame(
    n_cells  = nrow(d),
    n_orphan = sum(d$status == 1, na.rm = TRUE),
    n_b0     = sum(d$status == 2, na.rm = TRUE),
    rcri_p10 = q(d$rcri, 0.10), rcri_p25 = q(d$rcri, 0.25), rcri_med = q(d$rcri, 0.50),
    rcri_p75 = q(d$rcri, 0.75), rcri_p90 = q(d$rcri, 0.90))))
  out$area_km2 <- as.numeric(st_area(ap_ea)) / 1e6
  out$lat <- suppressWarnings(st_coordinates(st_point_on_surface(st_geometry(ap_sf)))[, 2])
  sf_out <- cbind(ap_ea[, COL_ID], out)
  f_shp <- file.path(DIR_OUT, paste0(tag, "_AP.shp"))
  write_sf(sf_out, f_shp, delete_dsn = file.exists(f_shp))
  write.csv(st_drop_geometry(sf_out), file.path(DIR_OUT, paste0(tag, "_AP.csv")), row.names = FALSE)
  invisible(sf_out)
}

tag_ra <- sprintf("RCRI_area_P%d", P)
res_a  <- rcri_ea(file.path(DIR_OUT, paste0(TAG_AA, "_EA_1km.tif")),
                  file.path(DIR_OUT, paste0(TAG_BA, "_EA_1km.tif")))
write_rcri(res_a, tag_ra)
ap_area <- resumen_ap(res_a, aps, tag_ra)


# 10 Complemento de replicación (cada AP suma 1) ----

if (RUN_REPLICATION) {

  # 10.1 mh_rep() por AP. Se salta si los Th_rep ya existen.
  n_chunks <- min(nrow(aps), NUM_CORES * 4)
  grp  <- split(seq_len(nrow(aps)), cut(seq_len(nrow(aps)), n_chunks, labels = FALSE))
  dirs <- file.path(DIR_RUN, sprintf("chunk_%03d", seq_along(grp)))
  g    <- rep(seq_along(grp), lengths(grp))
  expected <- file.path(dirs[g], "Analogues", paste0("Th_rep_", tolower(aps[[COL_ID]]), ".tif"))

  if (!all(file.exists(expected))) {
    cl <- makeCluster(NUM_CORES)
    registerDoParallel(cl)
    res_par <- foreach(i = seq_len(nrow(aps)), .packages = c("sf", "terra", "ClimaRep"),
                       .errorhandling = "pass") %dopar% {
                         if (file.exists(expected[i])) return("skip")
                         wtmp <- file.path(DIR_TMP, paste0("w", Sys.getpid()))
                         dir.create(wtmp, recursive = TRUE, showWarnings = FALSE)
                         terra::terraOptions(tempdir = wtmp, progress = 0, memfrac = 0.5 / NUM_CORES)
                         suppressMessages(ClimaRep::mh_rep(
                           polygon = aps[i, ], col_name = COL_ID, climate_variables = terra::rast(CLIM_VIF),
                           study_area = sa_box, th = TH, dir_output = dirs[g[i]], save_raw = FALSE))
                         invisible(gc())
                         unlink(list.files(wtmp, full.names = TRUE))
                         if (file.exists(expected[i])) "ok" else "sin salida"
                       }
    stopCluster(cl); registerDoSEQ()
    err <- vapply(res_par, inherits, logical(1), "error")
    if (any(err)) print(data.frame(id = aps[[COL_ID]][err],
                                   msg = vapply(res_par[err], conditionMessage, "")))
  }
  unlink(file.path(dirs, "Charts"), recursive = TRUE)
  ok_ap <- file.exists(expected)
  write.csv(data.frame(id = aps[[COL_ID]][!ok_ap]), file.path(DIR_OUT, "sites_without_output.csv"),
            row.names = FALSE)
  message("Sitios con salida: ", sum(ok_ap), " de ", length(ok_ap))

  # 10.2 Suma de Th_rep (equivale a rep_overlay(change = FALSE) si son 0/1)
  f_all <- expected[ok_ap]
  vals <- lapply(f_all, function(x) {
    fr <- try(terra::freq(terra::rast(x)), silent = TRUE)
    if (inherits(fr, "try-error")) NULL else fr$value
  })
  bad <- vapply(vals, is.null, logical(1))
  if (any(bad)) stop("Th_rep ilegibles: ", paste(basename(f_all[bad]), collapse = ", "))
  if (!all(unique(unlist(vals)) %in% c(0, 1))) stop("Los Th_rep no son binarios 0/1.")

  terraOptions(memfrac = 0.6)
  apcrs <- app(rast(f_all), "sum", na.rm = TRUE)
  # mh_rep() recorta el clima a sa_box, una celda más estrecho por cada lado
  clim1 <- rast(CLIM_VIF)[[1]]
  stopifnot(isTRUE(all.equal(res(apcrs), res(clim1))),
            max(abs(origin(apcrs) - origin(clim1))) < 1e-6)
  apcrs <- mask(extend(crop(apcrs, clim1, snap = "near"), clim1), clim1)
  stopifnot(compareGeom(apcrs, clim1),
            global(!is.na(clim1) & is.na(apcrs), "sum")[[1]] == 0)
  names(apcrs) <- "analogue_count"
  f_ap <- file.path(DIR_OUT, sprintf("APCRS_P%d_WGS84.tif", P))
  writeRaster(apcrs, f_ap, datatype = "INT4S", overwrite = TRUE, gdal = "COMPRESS=LZW")
  apcrs_ea <- mask(project(rast(f_ap), tpl, method = "near"), sa_ea)
  writeRaster(apcrs_ea, file.path(DIR_OUT, sprintf("APCRS_P%d_EA_1km.tif", P)),
              datatype = "INT4S", overwrite = TRUE, gdal = "COMPRESS=LZW")

  # 10.3 RCRI de replicación con el BCRS entero
  tag_rr <- sprintf("RCRI_rep_P%d", P)
  res_r  <- rcri_ea(file.path(DIR_OUT, sprintf("APCRS_P%d_EA_1km.tif", P)),
                    file.path(DIR_OUT, paste0(TAG_BI, "_EA_1km.tif")))
  write_rcri(res_r, tag_rr)
  ap_rep <- resumen_ap(res_r, aps[ok_ap, ], tag_rr)
}


# Se vacía la carpeta temporal sin borrarla, porque sigue siendo el tempdir de
# terra en la sesión.
unlink(list.files(DIR_TMP, full.names = TRUE), recursive = TRUE)
message("Terminado | ", round(as.numeric(difftime(Sys.time(), t0, units = "hours")), 2), " h")
