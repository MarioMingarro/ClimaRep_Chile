# Comparación entre la red actual y un escenario con anillos de amortiguación.
#
# Requiere haber ejecutado APCRS_BCRS_Chile.R con BUFFER_KM = 0 y con
# BUFFER_KM = ESC_KM (misma escala y umbral). Solo lee salidas.
#
# Medidas absolutas, comparables entre escenarios:
#   superficie protegida y añadida, huérfanos rellenados y nuevos, climas que
#   dejan de depender de una sola AP, factor de aumento del respaldo (APCRS en
#   km² de fragmentos) y eficiencia (km² de huérfanos rellenados por cada
#   1 000 km² anadidos).
# El RCRI es relativo (compara cuotas): si el anillo añade más superficie en
# unas zonas que en otras, el RCRI de las demás baja aunque su protección no
# cambie. Se incluye solo como referencia.
# Distancia de los huérfanos actuales a la AP más próxima: qué parte está al
# alcance de un anillo y qué parte exige áreas nuevas.
#
# El anillo no es superficie protegida: el escenario supone que las zonas de
# amortiguación se gestionan con criterios de conservación.

PATH_COMUN <- "Chile_comun.R"   # ruta al archivo común (Chile_config.R y Chile_unidades.R en la misma carpeta)
source(file.path(dirname(PATH_COMUN), "Chile_config.R"), encoding = "UTF-8")
source(file.path(dirname(PATH_COMUN), "Chile_unidades.R"), encoding = "UTF-8")
t0 <- Sys.time()

ESC_KM <- 5      # anchura del anillo del escenario que se compara
D0 <- dir_escenario(0)
D1 <- dir_escenario(ESC_KM)
DIR_CMP <- file.path(D1, "comparacion")
dir.create(DIR_CMP, recursive = TRUE, showWarnings = FALSE)
lab <- sprintf("B%gkm", ESC_KM)

capa <- function(d, f) {
  p <- file.path(d, f)
  if (!file.exists(p)) stop("Falta ", p, ". Ejecuta APCRS_BCRS_Chile.R en ese escenario.")
  rast(p)
}
f_st  <- paste0(TAG_RA, "_status_EA_1km.tif"); f_rc <- paste0(TAG_RA, "_EA_1km.tif")
f_red <- sprintf("redundancia_AP_%dkm_P%d_EA_1km.tif", UNIT_KM, P); f_a <- paste0(TAG_AA, "_EA_1km.tif")


# 1 Capas de los dos escenarios ----

st0 <- capa(D0, f_st); st1 <- capa(D1, f_st)
stopifnot(compareGeom(st0, st1))
cell <- which(!is.na(values(st0, mat = FALSE)))
stopifnot(identical(cell, which(!is.na(values(st1, mat = FALSE)))))
v <- function(r) values(r, mat = FALSE)[cell]
s0 <- v(st0); s1 <- v(st1)
r0 <- v(capa(D0, f_red)); r1 <- v(capa(D1, f_red))
a0 <- v(capa(D0, f_a));   a1 <- v(capa(D1, f_a))
k0 <- v(capa(D0, f_rc));  k1 <- v(capa(D1, f_rc))
h0 <- s0 == 1; h1 <- s1 == 1

aps0 <- read_sf(PATH_AP) |> st_transform(CRS_EA)
f_ap1 <- file.path(D1, "APs_escenario.gpkg")
if (!file.exists(f_ap1)) stop("Falta ", f_ap1, ".")
aps1 <- read_sf(f_ap1) |> st_transform(CRS_EA)
p0 <- v(rasterize(vect(aps0), st0, field = 1, background = 0, touches = FALSE)) == 1
p1 <- v(rasterize(vect(aps1), st0, field = 1, background = 0, touches = FALSE)) == 1

# Distancia de cada celda a la AP actual más próxima (km).
r_ap0 <- rasterize(vect(aps0), st0, field = 1, touches = FALSE)
dist0 <- v(distance(r_ap0)) / 1000

u_rm <- lee_regiones()
v_reg <- asigna_unidades(u_rm$reg, st0, cell)$r[cell]
v_mz  <- asigna_unidades(u_rm$mz,  st0, cell)$r[cell]
kp <- lee_koppen()
v_kp <- if (is.null(kp)) NULL else asigna_unidades(kp, st0, cell)$r[cell]

# Clases de cambio.
cambio <- ifelse(h0 & !h1, 1L,                     # huérfano rellenado
          ifelse(h0 & h1, 2L,                       # sigue huérfano
          ifelse(!h0 & h1, 3L,                      # nuevo huérfano
          ifelse(r0 == 1 & r1 >= 2, 4L, 5L))))      # deja de depender de una sola AP / sin cambio de clase
message(sprintf(paste("Celdas %d | protegidas %d -> %d | huérfanos %d -> %d (rellenados %d, nuevos %d)",
                      "| una sola AP %d -> %d"),
                length(cell), sum(p0), sum(p1), sum(h0), sum(h1), sum(cambio == 1), sum(cambio == 3),
                sum(r0 == 1), sum(r1 == 1)))


# 2 Tablas por unidad ----

resume_cmp <- function(g, nombres) {
  do.call(rbind, lapply(seq_along(nombres), function(i) {
    s <- which(g == i)
    anadido <- sum(p1[s]) - sum(p0[s])
    ok <- s[a0[s] > 0 & a1[s] > 0]
    data.frame(unidad = nombres[i], km2 = length(s),
               km2_protegidos_actual = sum(p0[s]), km2_protegidos_escenario = sum(p1[s]),
               km2_anadidos = anadido,
               pct_protegido_actual = 100 * mean(p0[s]), pct_protegido_escenario = 100 * mean(p1[s]),
               km2_huerfanos_actual = sum(h0[s]), km2_huerfanos_escenario = sum(h1[s]),
               km2_rellenados = sum(cambio[s] == 1), km2_nuevos_huerfanos = sum(cambio[s] == 3),
               pct_huerfanos_rellenados = if (sum(h0[s])) 100 * sum(cambio[s] == 1) / sum(h0[s]) else NA_real_,
               pct_una_AP_actual = 100 * mean(r0[s] == 1), pct_una_AP_escenario = 100 * mean(r1[s] == 1),
               factor_respaldo_mediano = if (length(ok)) 10^median(log10(a1[ok] / a0[ok])) else NA_real_,
               rellenados_por_1000km2 = if (anadido > 0) 1000 * sum(cambio[s] == 1) / anadido else NA_real_,
               pct_huerfanos_a_menos_del_anillo = if (sum(h0[s])) 100 * mean(dist0[s][h0[s]] <= ESC_KM) else NA_real_,
               RCRI_mediana_actual = if (any(s0[s] == 0)) median(k0[s][s0[s] == 0]) else NA_real_,
               RCRI_mediana_escenario = if (any(s1[s] == 0)) median(k1[s][s1[s] == 0]) else NA_real_)
  }))
}
orden_mz <- c("Norte", "Centro", "Centro Sur", "Sur", "Austral")
t_mz  <- resume_cmp(v_mz, u_rm$mz[[COL_MACROZONA]])
t_mz  <- rbind(t_mz[order(match(t_mz$unidad, orden_mz)), ], resume_cmp(rep(1L, length(cell)), "Chile"))
t_reg <- resume_cmp(v_reg, u_rm$reg[[COL_REGION]])
reg_mz <- vapply(seq_len(nrow(u_rm$reg)), function(i) {
  x <- v_mz[v_reg == i]; if (!length(x)) NA_integer_ else as.integer(names(which.max(table(x))))
}, integer(1))
t_reg <- cbind(t_reg[1], macrozona = u_rm$mz[[COL_MACROZONA]][reg_mz], t_reg[-1])
t_reg <- t_reg[order(match(t_reg$macrozona, orden_mz), -t_reg$km2), ]
names(t_mz)[1] <- "macrozona"; names(t_reg)[1] <- "region"
write.csv(t_mz,  file.path(DIR_CMP, sprintf("comparacion_macrozonas_%s.csv", lab)), row.names = FALSE)
write.csv(t_reg, file.path(DIR_CMP, sprintf("comparacion_regiones_%s.csv", lab)), row.names = FALSE)
if (!is.null(v_kp)) {
  t_kp <- resume_cmp(v_kp, kp$tipo)
  t_kp <- cbind(tipo = t_kp$unidad, descripcion = kp$descripcion, t_kp[, -1])
  t_kp <- t_kp[order(-t_kp$km2_huerfanos_actual, -t_kp$km2), ]
  write.csv(t_kp, file.path(DIR_CMP, sprintf("comparacion_koppen_%s.csv", lab)), row.names = FALSE)
}
op <- options(width = 220)
print(t_mz[, c("macrozona", "pct_protegido_actual", "pct_protegido_escenario", "km2_anadidos",
               "km2_huerfanos_actual", "km2_huerfanos_escenario", "pct_huerfanos_rellenados",
               "km2_nuevos_huerfanos", "pct_una_AP_actual", "pct_una_AP_escenario",
               "factor_respaldo_mediano", "rellenados_por_1000km2", "pct_huerfanos_a_menos_del_anillo")],
      digits = 3, row.names = FALSE)
options(op)


# 3 Distancia de los huérfanos actuales a la AP más próxima ----

cortes <- unique(c(0, ESC_KM, 10, 25, 50, 100, Inf))
etq <- paste0(head(cortes, -1), "–", tail(cortes, -1), " km"); etq[length(etq)] <- paste0("> ", cortes[length(cortes) - 1], " km")
banda <- cut(dist0[h0], cortes, labels = etq, include.lowest = TRUE)
t_dist <- data.frame(distancia = etq,
                     km2_huerfanos = as.vector(table(banda)),
                     km2_rellenados = as.vector(tapply(cambio[h0] == 1, banda, sum)))
t_dist$km2_rellenados[is.na(t_dist$km2_rellenados)] <- 0
t_dist$pct_del_total <- 100 * t_dist$km2_huerfanos / sum(h0)
write.csv(t_dist, file.path(DIR_CMP, sprintf("distancia_huerfanos_%s.csv", lab)), row.names = FALSE)
print(t_dist, digits = 3, row.names = FALSE)


# 4 Superficie de cada anillo por AP ----

f_an <- file.path(D1, sprintf("anillos_AP_B%gkm.csv", ESC_KM))
if (file.exists(f_an)) {
  t_an <- read.csv(f_an)
  t_an$pct_aumento <- 100 * t_an$area_anillo_km2 / t_an$area_nucleo_km2
  write.csv(t_an[order(-t_an$area_anillo_km2), ], file.path(DIR_CMP, sprintf("anillos_AP_%s.csv", lab)), row.names = FALSE)
}


# 5 Capa y figuras ----

etq_c <- c("Huérfano rellenado", "Sigue huérfano", "Nuevo huérfano", "Deja de depender de una sola AP", "Sin cambio de clase")
col_c <- c("#2166AC", "#B2182B", "#1F2A2E", "#E3B448", "#E6E3DA")
r_c <- rast(st0); vv <- rep(NA_real_, ncell(r_c)); vv[cell] <- cambio; values(r_c) <- vv
levels(r_c) <- data.frame(id = 1:5, cambio = etq_c)
coltab(r_c) <- data.frame(value = 1:5, col = col_c)
writeRaster(r_c, file.path(DIR_CMP, sprintf("cambio_%s_EA_1km.tif", lab)), datatype = "INT1U",
            overwrite = TRUE, gdal = "COMPRESS=LZW")

# Mapa del norte, centro y centro sur, donde están los huérfanos.
zona <- u_rm$mz[u_rm$mz[[COL_MACROZONA]] %in% c("Norte", "Centro", "Centro Sur"), ]
ez <- st_bbox(zona)
png(file.path(DIR_CMP, sprintf("mapa_cambio_%s.png", lab)), width = 1500, height = 2400, res = 250, bg = "#F6F4EE")
plot(crop(r_c, ext(ez["xmin"] - 20000, ez["xmax"] + 20000, ez["ymin"] - 10000, ez["ymax"] + 10000)),
     axes = FALSE, box = FALSE, mar = c(1, 1, 1, 9), plg = list(title = "Con el anillo", cex = 0.75))
plot(st_geometry(u_rm$reg), add = TRUE, border = "grey45", lwd = 0.4)
plot(st_geometry(zona), add = TRUE, border = "#1F2A2E", lwd = 1.2)
plot(st_geometry(aps0), add = TRUE, border = "#1F2A2E", lwd = 0.5)
dev.off()

# Huérfanos actuales y con el anillo por macrozona, y eficiencia.
tm <- t_mz[t_mz$macrozona %in% orden_mz, ]
png(file.path(DIR_CMP, sprintf("huerfanos_macrozonas_%s.png", lab)), width = 2400, height = 1300, res = 200, bg = "#F6F4EE")
par(mar = c(4.5, 9, 3, 12), bg = "#F6F4EE")
m <- rbind(tm$km2_huerfanos_escenario, tm$km2_rellenados)
b <- barplot(m[, nrow(tm):1] / 1000, horiz = TRUE, names.arg = rev(tm$macrozona), las = 1, border = NA,
             col = c("#B2182B", "#2166AC"), xlab = "Miles de km² de climas sin protección en la red actual",
             cex.names = 1.3, cex.axis = 1.2, cex.lab = 1.2, space = 0.35)
txt <- ifelse(is.na(tm$rellenados_por_1000km2), "", sprintf("%.0f km² por cada 1 000 anadidos", tm$rellenados_por_1000km2))
text(par("usr")[2], b, rev(txt), adj = 0, xpd = TRUE, cex = 1.05)
legend("topright", c("Siguen sin protección", sprintf("Rellenados por el anillo de %g km", ESC_KM)),
       fill = c("#B2182B", "#2166AC"), border = NA, bty = "n", cex = 1.1, inset = c(0, -0.12), xpd = TRUE, horiz = TRUE)
dev.off()

# Huérfanos actuales por distancia a la AP más próxima.
png(file.path(DIR_CMP, sprintf("distancia_huerfanos_%s.png", lab)), width = 2400, height = 1200, res = 200, bg = "#F6F4EE")
par(mar = c(5, 6, 3, 2), bg = "#F6F4EE")
md <- rbind(t_dist$km2_huerfanos - t_dist$km2_rellenados, t_dist$km2_rellenados)
barplot(md / 1000, names.arg = t_dist$distancia, col = c("#B2182B", "#2166AC"), border = NA,
        ylab = "Miles de km²", xlab = "Distancia a la área protegida más próxima",
        cex.names = 1.2, cex.axis = 1.2, cex.lab = 1.2)
legend("topright", c("Siguen sin protección", sprintf("Rellenados por el anillo de %g km", ESC_KM)),
       fill = c("#B2182B", "#2166AC"), border = NA, bty = "n", cex = 1.1)
dev.off()

message("Terminado | ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
