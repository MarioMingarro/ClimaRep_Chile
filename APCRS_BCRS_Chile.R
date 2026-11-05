# Representatividad climática de las AP de Chile por superficie.
# Fondo (BCRS) y AP (APCRS) con el mismo estimador. Cada cuadrícula del fondo es
# la tierra que contiene y cada fragmento protegido es la parte de una sola AP
# dentro de una cuadrícula. Cada unidad o fragmento tiene su envolvente P99
# (regla de celdas de ClimaRep) y suma su superficie / UNIT_KM^2 en las celdas
# de su envolvente. RCRI por celda en la rejilla de 1 km, con los huérfanos
# fuera del RCRI y en la capa de estado.
#
# Cada AP aporta la suma de las envolventes de sus fragmentos. De ahí salen,
# con el mismo cálculo:
#   APCRS        suma de los aportes de todas las AP.
#   Redundancia  número de AP distintas cuyo aporte incluye cada clima.
#   Por AP       RCRI observado y RCRI sin la propia AP (APCRS menos su aporte),
#                cuota propia del respaldo de sus climas y celdas exclusivas
#                (climas que solo respalda esa AP).
# En las cuadrículas con solapes, quitar una AP deja a las demás con su parte
# del reparto y no con toda la superficie solapada; la superficie solapada se
# registra en resumen_ejecucion.

PATH_COMUN <- "Chile_comun.R"   # ruta al archivo común
t0 <- Sys.time()
source(PATH_COMUN, encoding = "UTF-8")


# 9 BCRS y aporte de cada AP ----

cl <- abre_cluster()
bcrs_a <- suma_envolventes(idx_u[units_u], w_u[units_u], cl)
message("BCRS listo | ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
ap_frags <- split(frags, factor(frag_ap[frags], levels = seq_len(nrow(aps))))
tareas   <- lapply(ap_frags, function(f) list(p = idx_f[f], w = w_f[f]))
aporte   <- parLapplyLB(cl, tareas, trozo_env)
stopCluster(cl)
message("APCRS listo | ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")

apcrs_a <- numeric(length(z1s))
redund  <- integer(length(z1s))
for (r in aporte) {
  apcrs_a[r$i] <- apcrs_a[r$i] + r$x
  redund[r$i]  <- redund[r$i] + 1L
}

b_ea   <- bcrs_a[p_ea]
a_ea   <- apcrs_a[p_ea]
red_ea <- redund[p_ea]
escribe_ea(b_ea, TAG_BA)
escribe_ea(a_ea, TAG_AA)
escribe_ea(red_ea, sprintf("redundancia_AP_%dkm_P%d", UNIT_KM, P), dtype = "INT2U")
anota(celdas_clima_BCRS_0 = sum(bcrs_a == 0), celdas_clima_APCRS_0 = sum(apcrs_a == 0),
      celdas_EA_redundancia_0 = sum(red_ea == 0), celdas_EA_redundancia_1 = sum(red_ea == 1),
      celdas_EA_redundancia_2omas = sum(red_ea >= 2), redundancia_max = max(red_ea))
rm(bcrs_a, apcrs_a, redund)


# 10 RCRI por superficie ----

res <- rcri_vec(a_ea, b_ea)
escribe_ea(res$r, TAG_RA)
escribe_ea(res$st, paste0(TAG_RA, "_status"), dtype = "INT1U")
anota(celdas_EA = n_ea, celdas_validas = sum(res$st == 0),
      pct_RCRI_positivo = 100 * mean(res$r > 0, na.rm = TRUE),
      huerfanos = sum(res$st == 1), pct_huerfanos = 100 * mean(res$st == 1),
      BCRS_0_con_APCRS = sum(res$st == 2),
      RCRI_min = min(res$r, na.rm = TRUE), RCRI_mediana = median(res$r, na.rm = TRUE),
      RCRI_max = max(res$r, na.rm = TRUE))


# 11 Resumen por AP ----
# Celdas por su centro, sin fracciones. med_loo se calcula sobre las celdas que
# no quedan huérfanas sin la AP; hay que leerla junto a frac_exclusivas.

q   <- function(x, p) if (all(is.na(x))) NA_real_ else unname(quantile(x, p, na.rm = TRUE, type = 7))
Sb  <- sum(b_ea)
TOL <- 1e-9 * max(a_ea)
dv  <- numeric(length(z1s))
tab <- do.call(rbind, lapply(seq_len(nrow(aps)), function(k) {
  j_all <- cells_ap[[k]]; j <- j_all[!is.na(j_all)]
  r <- res$r[j]; s <- res$st[j]
  dv[] <- 0; dv[aporte[[k]]$i] <- aporte[[k]]$x
  a_k <- a_ea - dv[p_ea]
  a_k[a_k < TOL] <- 0
  ao <- a_ea[j]; ak <- a_k[j]; bj <- b_ea[j]
  r_loo <- ifelse(ak > 0 & bj > 0, RCRI_LOG((ak / sum(a_k)) / (bj / Sb)), NA_real_)
  cuota <- ifelse(ao > 0, (ao - ak) / ao, NA_real_)
  data.frame(n_cells = length(j_all), n_orphan = sum(s == 1L), n_b0 = sum(s == 2L),
             rcri_p10 = q(r, 0.10), rcri_p25 = q(r, 0.25), rcri_med = q(r, 0.50),
             rcri_p75 = q(r, 0.75), rcri_p90 = q(r, 0.90),
             med_loo = q(r_loo, 0.50), cuota_propia = q(cuota, 0.50),
             n_exclusivas = sum(ao > 0 & ak == 0), frac_exclusivas = mean(ao > 0 & ak == 0),
             redund_med = q(red_ea[j], 0.50),
             pct_masa_propia = 100 * sum(dv[p_ea]) / sum(a_ea),
             n_fragmentos = length(ap_frags[[k]]))
}))
tab$autocontrib <- tab$rcri_med - tab$med_loo
tab$area_km2 <- as.numeric(st_area(ap_ea2)) / 1e6
tab$lat      <- suppressWarnings(st_coordinates(st_point_on_surface(st_geometry(aps)))[, 2])
ap_out <- cbind(ap_ea2[, COL_ID], tab)
f_gpkg <- file.path(DIR_OUT, paste0(TAG_RA, "_AP.gpkg"))
write_sf(ap_out, f_gpkg, delete_dsn = file.exists(f_gpkg))
write.csv(st_drop_geometry(ap_out), file.path(DIR_OUT, paste0(TAG_RA, "_AP.csv")), row.names = FALSE)
anota(AP_mediana_positiva = sum(tab$rcri_med > 0, na.rm = TRUE),
      AP_med_loo_positiva = sum(tab$med_loo > 0, na.rm = TRUE),
      huerfanos_dentro_AP = sum(tab$n_orphan), AP_con_exclusivas = sum(tab$n_exclusivas > 0),
      celdas_exclusivas = sum(tab$n_exclusivas))


# 12 Resumen de la ejecución ----

anota(horas = as.numeric(difftime(Sys.time(), t0, units = "hours")))
rr <- data.frame(medida = names(RES_RUN), valor = vapply(RES_RUN, function(x) signif(x, 6), numeric(1)),
                 row.names = NULL)
write.csv(rr, file.path(DIR_OUT, sprintf("resumen_ejecucion_P%d_%dkm.csv", P, UNIT_KM)), row.names = FALSE)
op <- options(scipen = 999); print(rr, row.names = FALSE); options(op)

unlink(list.files(DIR_TMP, full.names = TRUE), recursive = TRUE)
message("Terminado | ", round(RES_RUN$horas, 2), " h")
