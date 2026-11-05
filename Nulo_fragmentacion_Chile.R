# Modelo nulo de fragmentación para el APCRS por superficie.
#
# Pregunta. Si la protección no seleccionara climas dentro de cada cuadrícula,
# ¿daría el estimador por fragmentos el mismo RCRI que la protección resumida a
# escala de cuadrícula? Se comparan tres APCRS frente al mismo BCRS_area.
#   OBS   fragmentos observados, igual que en APCRS_BCRS_Chile.R.
#   REF   cada fragmento aporta su peso w_f con la envolvente de su cuadrícula
#         completa. Es lo esperable sin selección dentro de la cuadrícula y sin
#         efecto de la fragmentación.
#   NULO  cada fragmento parcial se sustituye por un bloque compacto con su mismo
#         número de celdas, semilla al azar dentro de su cuadrícula y el mismo
#         peso. Los fragmentos con todas las celdas de su cuadrícula no cambian.
#         Los fragmentos son por AP, así que en una cuadrícula con varias AP
#         cada una se recoloca por separado.
# NULO − REF es el sesgo propio del estimador por fragmentos (lo que se quiere
# saber que es pequeño). OBS − NULO es la selección climática real dentro de las
# cuadrículas. REF por sí solo recoge la geografía de la red a la escala de la malla,
# por eso no se espera que el NULO quede centrado en cero.
#
# Limitaciones del nulo. Los bloques son compactos y de una sola parte, mientras
# que muchos fragmentos observados son multiparte. Con semilla uniforme y vecinos
# más próximos, las celdas del centro de la cuadrícula entran algo más a menudo
# que las de las esquinas.
#
# Usa la configuración de Chile_comun.R y depende de una ejecución de
# APCRS_BCRS_Chile.R con esa misma configuración (lee BCRS_area). La coherencia
# se comprueba recalculando el APCRS y el RCRI observados y comparándolos con
# los archivos.


PATH_COMUN <- "Chile_comun.R"   # ruta al archivo común
t0 <- Sys.time()
source(PATH_COMUN, encoding = "UTF-8")

N_REP       <- 30          # réplicas del nulo
SEED        <- 20261104    # semilla base, la réplica r usa SEED + r
RECALC_BCRS <- FALSE       # TRUE recalcula el BCRS_area en lugar de leerlo

DIR_NULL <- file.path(DIR_OUT, "nulo_fragmentacion")
dir.create(DIR_NULL, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(file.path(DIR_OUT, paste0(TAG_BA, "_EA_1km.tif"))) && !RECALC_BCRS)
  stop("Falta ", TAG_BA, " en ", DIR_OUT, ". Ejecuta antes APCRS_BCRS_Chile.R con UNIT_KM = ", UNIT_KM, ".")

n_f    <- lengths(idx_f[frags])
n_u    <- lengths(idx_u[frag_sq[frags]])
full_f <- n_f == n_u
message(sprintf("Fragmentos %d | con todas las celdas de su cuadrícula %d | parciales %d",
                length(frags), sum(full_f), sum(!full_f)))


# 8 BCRS por superficie ----
# Un solo grupo de nodos para el BCRS (si se recalcula), OBS, REF y las réplicas.

cl <- abre_cluster()
if (RECALC_BCRS) {
  b_ea <- suma_envolventes(idx_u[units_u], w_u[units_u], cl)[p_ea]
} else {
  b_ea <- lee_ea(file.path(DIR_OUT, paste0(TAG_BA, "_EA_1km.tif")))
}


# 9 APCRS observado y de referencia, y comprobación de coherencia ----
# REF da a cada fragmento la envolvente de su cuadrícula completa. base_full es
# el aporte fijo de los fragmentos con todas las celdas de su cuadrícula.

a_obs     <- suma_envolventes(idx_f[frags], w_f[frags], cl)
a_ref     <- suma_envolventes(idx_u[frag_sq[frags]], w_f[frags], cl)
base_full <- suma_envolventes(idx_u[frag_sq[frags[full_f]]], w_f[frags[full_f]], cl)
a_obs_ea <- a_obs[p_ea]
a_ref_ea <- a_ref[p_ea]
rm(a_obs, a_ref)

a_file <- lee_ea(file.path(DIR_OUT, paste0(TAG_AA, "_EA_1km.tif")))
d_a <- max(abs(a_obs_ea - a_file) / pmax(1, abs(a_file)))
obs <- rcri_vec(a_obs_ea, b_ea)
r_file <- values(rast(file.path(DIR_OUT, paste0(TAG_RA, "_EA_1km.tif"))), mat = FALSE)[cell_ea]
d_r <- max(abs(obs$r - r_file), na.rm = TRUE)
same_na <- identical(is.na(obs$r), is.na(r_file))
message(sprintf("Coherencia con la ejecución principal | APCRS dif. rel. máx %.2e | RCRI dif. máx %.2e | NA iguales %s",
                d_a, d_r, same_na))
if (d_a > 1e-5 || d_r > 1e-4 || !same_na)
  stop("El APCRS o el RCRI recalculados no coinciden con los archivos de DIR_OUT.")
rm(a_file, r_file)
ref <- rcri_vec(a_ref_ea, b_ea)


# 10 Réplicas del nulo ----
# Bloque compacto, semilla uniforme entre las celdas de la cuadrícula y sus n
# vecinas más próximas en LAEA (empates al azar). Cada réplica se guarda en disco
# y se reutiliza si su huella coincide con la configuración actual.

part <- frags[!full_f]
FR <- lapply(part, function(u) {
  ci <- cells_sq[[frag_sq[[u]]]]
  list(ci = ci, xy = xy_ea[ci, , drop = FALSE], n = length(idx_f[[u]]), w = w_f[[u]])
})
FP <- c(n_ok = length(ok), n_ea = n_ea, n_frag = length(frags), n_part = length(part),
        sum_w = sum(w_f[frags]), TH = TH, UNIT_KM = UNIT_KM, SEED = SEED)

null_rep <- function(r) {
  set.seed(SEED + r)
  v <- base_full
  for (f in FR) {
    m <- nrow(f$xy)
    s <- sample.int(m, 1L)
    d <- (f$xy[, 1] - f$xy[s, 1])^2 + (f$xy[, 2] - f$xy[s, 2])^2
    b <- order(d, runif(m))[seq_len(f$n)]
    h <- env_hits(pos[f$ci[b]])$hit
    v[h] <- v[h] + f$w
  }
  saveRDS(list(fp = FP, a_ea = v[p_ea]), file.path(DIR_NULL, sprintf("rep_%03d.rds", r)))
  r
}

f_rep  <- file.path(DIR_NULL, sprintf("rep_%03d.rds", seq_len(N_REP)))
valida <- vapply(f_rep, function(f) file.exists(f) && isTRUE(all.equal(readRDS(f)$fp, FP)), logical(1))
falta  <- which(!valida)
if (length(falta)) {
  message("Réplicas por calcular: ", length(falta))
  t1 <- Sys.time()
  clusterExport(cl, c("pos", "FR", "base_full", "p_ea", "SEED", "FP", "DIR_NULL"), envir = globalenv())
  invisible(parLapply(cl, falta, null_rep))
  message("Réplicas en ", round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1), " min")
}
stopCluster(cl)


# 11 RCRI del nulo por celda y por AP ----

ap_ea <- ap_ea2

R_null <- matrix(NA_real_, n_ea, N_REP)
orf    <- integer(n_ea)
M_null <- matrix(NA_real_, nrow(ap_ea), N_REP)
n_orf  <- integer(N_REP)
n_val  <- integer(N_REP)
for (r in seq_len(N_REP)) {
  x <- rcri_vec(readRDS(f_rep[r])$a_ea, b_ea)
  R_null[, r] <- x$r
  orf <- orf + (x$st == 1L)
  n_orf[r] <- sum(x$st == 1L)
  n_val[r] <- sum(x$st == 0L)
  M_null[, r] <- med_ap(x$r)
}
n_ok_r  <- rowSums(!is.na(R_null))
nu_mean <- rowSums(R_null, na.rm = TRUE) / n_ok_r
nu_sd   <- sqrt(pmax(rowSums(R_null^2, na.rm = TRUE) - n_ok_r * nu_mean^2, 0) / pmax(n_ok_r - 1, 1))
nu_mean[n_ok_r == 0] <- NA_real_
nu_sd[n_ok_r < 2]    <- NA_real_
orf_frec <- orf / N_REP
rm(R_null); invisible(gc())

med_obs <- med_ap(obs$r)
f_csv <- file.path(DIR_OUT, paste0(TAG_RA, "_AP.csv"))
if (file.exists(f_csv)) {
  cv <- read.csv(f_csv)
  stopifnot(identical(as.character(cv[[COL_ID]]), as.character(ap_ea[[COL_ID]])))
  message(sprintf("Medianas por AP frente al CSV principal, dif. máx %.2e",
                  max(abs(cv$rcri_med - med_obs), na.rm = TRUE)))
}


# 12 Salidas ----

escribe_ea(ref$r, sprintf("RCRI_ref_P%d", P), DIR_NULL)
escribe_ea(ref$st, sprintf("RCRI_ref_P%d_status", P), DIR_NULL, "INT1U")
escribe_ea(nu_mean, sprintf("RCRI_nulo_media_P%d", P), DIR_NULL)
escribe_ea(nu_sd, sprintf("RCRI_nulo_sd_P%d", P), DIR_NULL)
escribe_ea(nu_mean - ref$r, sprintf("sesgo_fragmentacion_P%d", P), DIR_NULL)     # NULO − REF
escribe_ea(obs$r - nu_mean, sprintf("seleccion_intracuadricula_P%d", P), DIR_NULL) # OBS − NULO
escribe_ea(orf_frec, sprintf("huerfanos_nulo_frec_P%d", P), DIR_NULL)

mn  <- rowMeans(M_null, na.rm = TRUE)
med_ref <- med_ap(ref$r)
dev <- abs(M_null - mn)
p_b <- (1 + rowSums(dev >= abs(med_obs - mn), na.rm = TRUE)) / (1 + rowSums(!is.na(M_null)))
tab <- data.frame(
  SITECODE   = ap_ea[[COL_ID]],
  area_km2   = as.numeric(st_area(ap_ea)) / 1e6,
  lat        = suppressWarnings(st_coordinates(st_point_on_surface(st_geometry(aps)))[, 2]),
  n_cells    = lengths(cells_ap),
  med_obs    = med_obs,
  med_ref    = med_ref,
  med_nulo   = mn,
  nulo_q025  = apply(M_null, 1, quantile, 0.025, na.rm = TRUE),
  nulo_q975  = apply(M_null, 1, quantile, 0.975, na.rm = TRUE),
  sesgo      = mn - med_ref,       # NULO − REF
  seleccion  = med_obs - mn,       # OBS − NULO
  p_bilateral = p_b)
names(tab)[1] <- COL_ID
write.csv(tab, file.path(DIR_NULL, sprintf("nulo_fragmentacion_AP_P%d.csv", P)), row.names = FALSE)


# 13 Resumen ----
# El diseño queda defendido si NULO − REF es pequeño frente a OBS − REF y si la
# relación de la mediana por AP con el área y la latitud aparece ya en REF.

rho <- function(x, y) suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))
v   <- !is.na(obs$r) & !is.na(ref$r) & !is.na(nu_mean)
res <- data.frame(
  medida = c("celdas válidas", "huérfanos", "celdas con RCRI > 0 (%)", "mediana RCRI",
             "rho mediana AP con área", "rho mediana AP con latitud"),
  OBS  = c(sum(obs$st == 0), sum(obs$st == 1), 100 * mean(obs$r > 0, na.rm = TRUE),
           median(obs$r, na.rm = TRUE), rho(med_obs, tab$area_km2), rho(med_obs, tab$lat)),
  REF  = c(sum(ref$st == 0), sum(ref$st == 1), 100 * mean(ref$r > 0, na.rm = TRUE),
           median(ref$r, na.rm = TRUE), rho(med_ref, tab$area_km2), rho(med_ref, tab$lat)),
  NULO = c(mean(n_val), mean(n_orf), 100 * mean(nu_mean > 0, na.rm = TRUE),
           median(nu_mean, na.rm = TRUE), rho(mn, tab$area_km2), rho(mn, tab$lat)))
print(res, digits = 3)

dif <- data.frame(
  comparacion = c("NULO − REF (sesgo)", "OBS − NULO (selección)", "OBS − REF (total)"),
  media = c(mean(nu_mean[v] - ref$r[v]), mean(obs$r[v] - nu_mean[v]), mean(obs$r[v] - ref$r[v])),
  sd    = c(sd(nu_mean[v] - ref$r[v]), sd(obs$r[v] - nu_mean[v]), sd(obs$r[v] - ref$r[v])),
  p05   = c(quantile(nu_mean[v] - ref$r[v], 0.05), quantile(obs$r[v] - nu_mean[v], 0.05),
            quantile(obs$r[v] - ref$r[v], 0.05)),
  p95   = c(quantile(nu_mean[v] - ref$r[v], 0.95), quantile(obs$r[v] - nu_mean[v], 0.95),
            quantile(obs$r[v] - ref$r[v], 0.95)),
  rho_celdas = c(rho(nu_mean[v], ref$r[v]), rho(obs$r[v], nu_mean[v]), rho(obs$r[v], ref$r[v])),
  sd_AP = c(sd(tab$sesgo, na.rm = TRUE), sd(tab$seleccion, na.rm = TRUE),
            sd(tab$med_obs - tab$med_ref, na.rm = TRUE)),
  rho_AP_area = c(rho(tab$sesgo, tab$area_km2), rho(tab$seleccion, tab$area_km2),
                  rho(tab$med_obs - tab$med_ref, tab$area_km2)))
print(dif, digits = 3)
message(sprintf("Huérfanos del nulo en todas las réplicas %d | en la mitad o más %d | huérfanos REF %d | OBS %d",
                sum(orf_frec == 1), sum(orf_frec >= 0.5), sum(ref$st == 1), sum(obs$st == 1)))
message(sprintf("AP con OBS fuera del intervalo 95 %% del nulo %d de %d",
                sum(tab$med_obs < tab$nulo_q025 | tab$med_obs > tab$nulo_q975, na.rm = TRUE), nrow(tab)))
write.csv(res, file.path(DIR_NULL, sprintf("nulo_resumen_P%d.csv", P)), row.names = FALSE)
write.csv(dif, file.path(DIR_NULL, sprintf("nulo_diferencias_P%d.csv", P)), row.names = FALSE)

unlink(list.files(DIR_TMP, full.names = TRUE), recursive = TRUE)
message("Terminado | ", round(as.numeric(difftime(Sys.time(), t0, units = "hours")), 2), " h")
