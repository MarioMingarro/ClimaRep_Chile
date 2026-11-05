# Preparación común a APCRS_BCRS_Chile.R y Nulo_fragmentacion_Chile.R, que lo
# cargan con source(). La configuración está en Chile_config.R, que también
# usan los scripts que solo leen salidas, así que nada puede desincronizarse.
# Cada escala (UNIT_KM) y umbral tiene su propia carpeta de salida; el clima es
# común.
# Las cifras de control de la ejecución se acumulan en RES_RUN y el script
# principal las guarda en resumen_ejecucion.csv.
#
# Fragmentos por AP y cuadrícula. Una envolvente nunca reúne AP distintas: cada
# fragmento es la parte de una sola AP dentro de una cuadrícula. La superficie
# solapada entre AP se reparte a partes iguales entre ellas y las celdas
# solapadas entran en la envolvente de cada una.


# 0 Configuración ----
# Vive en Chile_config.R, en la misma carpeta que este archivo.

if (!exists("PATH_COMUN")) PATH_COMUN <- "Chile_comun.R"
source(file.path(dirname(PATH_COMUN), "Chile_config.R"), encoding = "UTF-8")


# 1 Área de estudio y rejilla equivalente ----

sa <- read_sf(PATH_STUDY_AREA) |> st_make_valid() |> st_transform(4326) |> st_union() |> st_make_valid()
sa <- st_sf(ID = 1L, geometry = sa)
sa_ea_sf <- st_make_valid(st_transform(sa, CRS_EA))
sa_ea    <- vect(sa_ea_sf)
e <- as.vector(ext(sa_ea))
tpl <- rast(xmin = floor(e[1] / RES_EA) * RES_EA, xmax = ceiling(e[2] / RES_EA) * RES_EA,
            ymin = floor(e[3] / RES_EA) * RES_EA, ymax = ceiling(e[4] / RES_EA) * RES_EA,
            resolution = RES_EA, crs = CRS_EA)


# 2 Clima ----
# Se construye una sola vez y lo comparten todas las escalas. BIO1, BIO2, BIO4,
# BIO12 y BIO15 sin transformar. El VIF es solo diagnóstico y se calcula si
# ClimaRep está instalado.

VARS <- c("bio1", "bio2", "bio4", "bio12", "bio15")
if (!file.exists(CLIM_VIF)) {
  f <- list.files(DIR_CHELSA, "bio[0-9]+.*\\.tif$", full.names = TRUE, ignore.case = TRUE)
  v <- sub("^bio0+([1-9])", "bio\\1", sub("^.*?(bio[0-9]+).*$", "\\1", tolower(basename(f)), perl = TRUE))
  clim <- rast(f[match(VARS, v)]); names(clim) <- VARS
  crs(clim) <- "EPSG:4326"
  sa_v <- vect(sa)
  clim <- mask(crop(clim, sa_v, snap = "out"), sa_v, touches = FALSE)
  if (requireNamespace("ClimaRep", quietly = TRUE)) {
    vif <- ClimaRep::vif_filter(clim, th = VIF_TH)
    message("VIF retendría: ", paste(vif$summary$kept_layers, collapse = ", "))
  }
  writeRaster(clim, CLIM_VIF, overwrite = TRUE,
              gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "BIGTIFF=IF_SAFER"))
}
clim <- rast(CLIM_VIF)
stopifnot(identical(names(clim), VARS))


# 3 Áreas protegidas ----

aps <- read_sf(PATH_AP_RAW) |> st_make_valid() |> st_transform(CRS_EA) |> st_make_valid()
aps <- aggregate(aps[, COL_ID], by = list(id = aps[[COL_ID]]), FUN = function(x) x[1])[, COL_ID]
aps <- suppressWarnings(st_intersection(aps, st_geometry(sa_ea_sf)))
aps <- st_collection_extract(aps, "POLYGON", warn = FALSE)
aps <- aggregate(aps[, COL_ID], by = list(id = aps[[COL_ID]]), FUN = function(x) x[1])[, COL_ID]
aps <- aps[as.numeric(st_area(aps)) > MIN_AREA_M2, ]
if (anyDuplicated(tolower(aps[[COL_ID]]))) stop("Hay ", COL_ID, " duplicados.")
if (!file.exists(PATH_AP)) write_sf(st_transform(aps, 4326), PATH_AP)

# Escenario con anillos. A cada AP se le añade la parte de su buffer que no es
# AP (ninguna, ni la propia ni otra) y que cae en el área de estudio. Los anillos
# de AP vecinas pueden solaparse entre sí; esa superficie se reparte después a
# partes iguales, como cualquier solape (bloque 6).
if (BUFFER_KM > 0) {
  ids    <- aps[[COL_ID]]
  nucleo <- st_union(st_geometry(aps))
  anillo <- suppressWarnings(st_buffer(aps[, COL_ID], BUFFER_KM * 1000))
  anillo <- suppressWarnings(st_difference(anillo, nucleo))
  anillo <- suppressWarnings(st_intersection(anillo, st_geometry(sa_ea_sf)))
  anillo <- st_collection_extract(st_make_valid(anillo), "POLYGON", warn = FALSE)
  t_anillo <- data.frame(id = ids, area_nucleo_km2 = as.numeric(st_area(aps)) / 1e6,
                         area_anillo_km2 = 0)
  a_an <- tapply(as.numeric(st_area(anillo)), anillo[[COL_ID]], sum) / 1e6
  t_anillo$area_anillo_km2[match(names(a_an), ids)] <- a_an
  names(t_anillo)[1] <- COL_ID
  dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)
  write.csv(t_anillo, file.path(DIR_OUT, sprintf("anillos_AP_B%gkm.csv", BUFFER_KM)), row.names = FALSE)
  aps <- rbind(aps[, COL_ID], anillo[, COL_ID])
  aps <- aggregate(aps[, COL_ID], by = list(id = aps[[COL_ID]]), FUN = function(x) x[1])[, COL_ID]
  aps <- st_make_valid(aps)
  aps <- aps[match(ids, aps[[COL_ID]]), ]
  stopifnot(identical(aps[[COL_ID]], ids))
  anota(area_nucleo_km2 = sum(t_anillo$area_nucleo_km2), area_anillos_km2 = sum(t_anillo$area_anillo_km2))
  message(sprintf("Anillos de %g km | núcleos %.0f km2 | anillos %.0f km2 (sin contar solapes entre anillos)",
                  BUFFER_KM, sum(t_anillo$area_nucleo_km2), sum(t_anillo$area_anillo_km2)))
  write_sf(st_transform(aps, 4326), PATH_AP_ESC, delete_dsn = file.exists(PATH_AP_ESC))
}
aps_ea <- aps
aps    <- st_transform(aps, 4326)
anota(n_AP = nrow(aps), area_AP_km2 = sum(as.numeric(st_area(aps_ea))) / 1e6)
message("Áreas protegidas > 5 km2: ", nrow(aps))


# 4 Clima blanqueado y asignación de celdas ----
# Cada celda pertenece a la cuadrícula que contiene su centro y está protegida
# si su centro cae dentro de alguna AP (regla de mask(crop()) de ClimaRep).
# D2 en espacio blanqueado con la covarianza de toda el área de estudio.

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
idx_u    <- split(pos, uid)                 # cuadrículas del fondo, orden de Zs
cells_sq <- split(seq_along(ok), uid)       # celdas de cada cuadrícula, orden de ok
anota(celdas_clima = length(ok), celdas_protegidas = sum(prot),
      cuadriculas = length(idx_u), cuadriculas_con_proteccion = length(unique(uid[prot])))
message(sprintf("Celdas con clima %d | protegidas %d | cuadrículas %d | con celdas protegidas %d",
                length(ok), sum(prot), length(idx_u), RES_RUN$cuadriculas_con_proteccion))


# 5 Índice de celdas y celdas de cada AP ----
# Índice de las celdas de clima en WGS84 (para pasar a la rejilla de 1 km y
# para extraer las celdas de cada AP). Celdas por su centro, como en
# rasterize(touches = FALSE). Una celda puede estar en varias AP si se solapan.

f_idx <- file.path(DIR_TMP, "indice_celdas_WGS84.tif")
idx_r <- rast(clim, nlyrs = 1)
vals  <- rep(NA_real_, ncell(idx_r)); vals[ok] <- seq_along(ok); values(idx_r) <- vals
writeRaster(idx_r, f_idx, datatype = "INT4S", overwrite = TRUE, gdal = "COMPRESS=LZW")
rm(idx_r, vals)

ex_w    <- terra::extract(rast(f_idx), vect(aps), ID = TRUE)
cells_k <- lapply(split(ex_w[[2]], factor(ex_w$ID, levels = seq_len(nrow(aps)))),
                  function(x) unique(x[!is.na(x)]))
rm(ex_w)
cnt   <- tabulate(unlist(cells_k), nbins = length(ok))
n_dif <- sum(prot != (cnt > 0))
if (n_dif > 0.001 * sum(prot)) stop("La extracción de celdas por AP no coincide con la rasterización.")
anota(celdas_dif_rasterizacion_extraccion = n_dif, celdas_en_varias_AP = sum(cnt > 1))


# 6 Malla, pesos y fragmentos por AP ----
# Pesos a ambos lados con la misma definición, el área en LAEA dividida entre
# UNIT_KM^2. Cada fragmento es la parte de una AP en una cuadrícula (pid =
# cuadrícula|AP). En las cuadrículas con varias AP, la superficie solapada se
# reparte a partes iguales.

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

g_sa <- st_geometry(sa_ea_sf)
hit  <- lengths(st_intersects(sq_u, g_sa)) > 0
full <- hit & lengths(st_within(sq_u, g_sa)) > 0
land <- setNames(numeric(nrow(sq_u)), units_u)
land[full] <- s_m^2
brd <- which(hit & !full)
if (length(brd)) {
  xl <- suppressWarnings(st_intersection(sq_u[brd, ], g_sa))
  al <- tapply(as.numeric(st_area(xl)), xl$uid, sum)
  land[names(al)] <- al
}
w_u <- pmin(land / s_m^2, 1)
rm(hit, full, land, brd)

sq_ap <- sq_u[lengths(st_intersects(sq_u, aps_ea)) > 0, "uid"]
xk <- suppressWarnings(st_intersection(sq_ap, aps_ea[, COL_ID]))
xk <- st_collection_extract(st_make_valid(xk), "POLYGON", warn = FALSE)
xk$k   <- match(xk[[COL_ID]], aps_ea[[COL_ID]])
xk$pid <- paste(xk$uid, xk$k, sep = "|")
xk$a   <- as.numeric(st_area(xk))
xk <- xk[xk$a > 0, ]
a_frag  <- tapply(xk$a, xk$pid, sum)
a_bruta <- sum(a_frag)
n_parts <- table(suppressWarnings(st_cast(st_cast(xk[, "pid"], "MULTIPOLYGON"), "POLYGON"))$pid)

# Reparto de los solapes en las cuadrículas con varias AP.
n_ap_sq  <- tapply(xk$k, xk$uid, function(z) length(unique(z)))
multi_sq <- names(n_ap_sq)[n_ap_sq > 1]
solape   <- 0
for (u in multi_sq) {
  s  <- xk[xk$uid == u, "k"]
  s  <- aggregate(s, by = list(kk = s$k), FUN = function(x) x[1])
  ov <- suppressWarnings(st_collection_extract(st_intersection(s[, "kk"]), "POLYGON", warn = FALSE))
  ar <- as.numeric(st_area(ov)); no <- ov$n.overlaps
  if (any(no > 1)) {
    solape <- solape + sum(ar[no > 1] * (no[no > 1] - 1))
    for (r in seq_len(nrow(s))) {
      has <- vapply(ov$origins, function(og) r %in% og, logical(1))
      a_frag[[paste(u, s$kk[r], sep = "|")]] <- sum(ar[has] / no[has])
    }
  }
}
w_f     <- pmin(a_frag / s_m^2, 1)
frag_a  <- a_frag

# Celdas de cada fragmento (orden de Zs).
idx_f <- unlist(lapply(seq_along(cells_k), function(k) {
  ck <- cells_k[[k]]
  if (!length(ck)) return(list())
  sp <- split(pos[ck], uid[ck])
  names(sp) <- paste(names(sp), k, sep = "|")
  sp
}), recursive = FALSE)

sin_geom  <- setdiff(names(idx_f), names(w_f))
sin_celda <- setdiff(names(w_f), names(idx_f))
frags     <- intersect(names(idx_f), names(w_f))
frag_sq   <- setNames(sub("\\|.*$", "", frags), frags)
frag_ap   <- setNames(as.integer(sub("^.*\\|", "", frags)), frags)
anota(cuadriculas_sin_tierra = sum(w_u == 0),
      fragmentos_celdas_sin_geometria = length(sin_geom),
      celdas_descartadas_sin_geometria = sum(lengths(idx_f[sin_geom])),
      fragmentos = length(frags), fragmentos_multiparte = sum(n_parts[frags] > 1),
      fragmentos_sin_celda = length(sin_celda), area_fragmentos_sin_celda_km2 = sum(frag_a[sin_celda]) / 1e6,
      cuadriculas_con_AP = length(unique(xk$uid)), cuadriculas_con_varias_AP = length(multi_sq),
      area_en_cuadriculas_con_varias_AP_km2 = sum(a_frag[sub("\\|.*$", "", names(a_frag)) %in% multi_sq]) / 1e6,
      superficie_solapada_km2 = solape / 1e6,
      area_fragmentos_km2 = sum(frag_a) / 1e6, fraccion_protegida = sum(w_f[frags]) / sum(w_u))
message(sprintf(paste("Fragmentos por AP %d (multiparte %d) | sin centro de celda %d (%.2f km2)",
                      "| con celdas y sin geometría %d (%d celdas) | cuadrículas con AP %d, con varias AP %d (%.0f km2)",
                      "| superficie solapada %.2f km2 | superficie de fragmentos %.0f km2 | suma de áreas de AP %.0f km2",
                      "| fracción protegida %.4f"),
                RES_RUN$fragmentos, RES_RUN$fragmentos_multiparte, RES_RUN$fragmentos_sin_celda,
                RES_RUN$area_fragmentos_sin_celda_km2, RES_RUN$fragmentos_celdas_sin_geometria,
                RES_RUN$celdas_descartadas_sin_geometria, RES_RUN$cuadriculas_con_AP,
                RES_RUN$cuadriculas_con_varias_AP, RES_RUN$area_en_cuadriculas_con_varias_AP_km2,
                RES_RUN$superficie_solapada_km2, RES_RUN$area_fragmentos_km2, RES_RUN$area_AP_km2,
                RES_RUN$fraccion_protegida))


# 7 Envolventes ----
# Regla de celdas de ClimaRep, validada contra mh_rep() en 60 cuadrículas y 55
# fragmentos sin ninguna diferencia. D2 exacto con poda por el primer eje.

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

# Suma de pesos sobre las envolventes de una lista de unidades, en paralelo.
# trozo_env() se define aquí, en el entorno global, para que al enviarla a los
# nodos no arrastre objetos grandes.
trozo_env <- function(t) {
  v <- numeric(length(z1s))
  for (k in seq_along(t$p)) {
    h <- env_hits(t$p[[k]])$hit
    v[h] <- v[h] + t$w[k]
  }
  i <- which(v != 0)
  list(i = i, x = v[i])
}
abre_cluster <- function(extra = character(0)) {
  cl <- makeCluster(NUM_CORES)
  clusterExport(cl, c("Zs", "z1s", "z2s", "TH", "env_hits", "trozo_env", extra), envir = globalenv())
  cl
}
suma_envolventes <- function(p_list, w, cl = NULL) {
  w <- as.numeric(w)
  n <- length(p_list)
  if (is.null(cl) || n < 4 * NUM_CORES) {
    res <- list(trozo_env(list(p = p_list, w = w)))
  } else {
    ch  <- split(seq_len(n), cut(seq_len(n), 4 * NUM_CORES, labels = FALSE))
    res <- parLapply(cl, lapply(ch, function(ix) list(p = p_list[ix], w = w[ix])), trozo_env)
  }
  v <- numeric(length(z1s))
  for (r in res) v[r$i] <- v[r$i] + r$x
  v
}


# 8 Rejilla equivalente y resúmenes por AP ----
# Correspondencia directa entre la rejilla de 1 km y las celdas de clima, por
# vecino más próximo y con la máscara del área de estudio, igual que
# project(..., method = "near") seguido de mask(). Un vector v en el orden de
# Zs pasa a la rejilla con v[p_ea].

idx_ea  <- mask(project(rast(f_idx), tpl, method = "near"), sa_ea)
v_ea    <- values(idx_ea, mat = FALSE)
cell_ea <- which(!is.na(v_ea))
i_ea    <- v_ea[cell_ea]
stopifnot(all(i_ea == round(i_ea)), min(i_ea) >= 1, max(i_ea) <= length(ok))
p_ea    <- pos[i_ea]
n_ea    <- length(cell_ea)
rm(v_ea)

escribe_ea <- function(x, nombre, dir = DIR_OUT, dtype = "FLT4S") {
  r <- rast(tpl); v <- rep(NA_real_, ncell(r)); v[cell_ea] <- x; values(r) <- v
  names(r) <- nombre
  gd <- if (dtype == "FLT4S") c("COMPRESS=DEFLATE", "PREDICTOR=3") else "COMPRESS=LZW"
  writeRaster(r, file.path(dir, paste0(nombre, "_EA_1km.tif")), datatype = dtype,
              overwrite = TRUE, gdal = gd)
}
lee_ea <- function(f) {
  r <- rast(f)
  stopifnot(compareGeom(r, tpl))
  v <- values(r, mat = FALSE)
  if (!identical(which(!is.na(v)), cell_ea)) stop("La máscara de ", basename(f), " no coincide.")
  v[cell_ea]
}

# RCRI por celda. Estado 0 válido, 1 huérfano (APCRS = 0), 2 con BCRS = 0.
rcri_vec <- function(a, b) {
  st <- ifelse(a == 0, 1L, ifelse(b == 0, 2L, 0L))
  r  <- RCRI_LOG((a / sum(a)) / (b / sum(b)))
  r[st != 0] <- NA_real_
  list(r = r, st = st)
}

# Celdas de cada AP en la rejilla de 1 km, por su centro y sin fracciones.
ap_ea2   <- st_transform(aps, CRS_EA)
ex       <- terra::extract(idx_ea, vect(ap_ea2), cells = TRUE, ID = TRUE)
cells_ap <- split(match(ex$cell, cell_ea), factor(ex$ID, levels = seq_len(nrow(ap_ea2))))
rm(ex)
med_ap <- function(r) vapply(cells_ap, function(j) {
  x <- r[j]; if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
}, numeric(1))
