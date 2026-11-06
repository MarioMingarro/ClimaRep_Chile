# Huérfanos, redundancia y RCRI por región y macrozona, y caracterización
# climática de los huérfanos.
#
# Lee las salidas de APCRS_BCRS_Chile.R (método único) sin recalcular nada.
#   Tablas por región y por macrozona con superficie, superficie protegida, AP,
#   huérfanos, climas respaldados por una sola AP, RCRI mediano y fracción de
#   territorio muy infrarrepresentado (RCRI < −1, menos de una décima parte
#   de lo esperado).
#   Tabla climática con los percentiles 10, 50 y 90 de cada variable en todo
#   Chile, en los climas con análogo y en los huérfanos, también por macrozona.
#   Capa de clases de redundancia (huérfano, 1 AP, 2 a 4 AP, 5 o más AP).
#   Tabla y gráfico por tipo de clima de Köppen.
#   Figuras (mapa de clases y espacio climático de los huérfanos coloreado por
#   tipo de Köppen).
# Cada celda de la rejilla de 1 km en LAEA es 1 km².

PATH_COMUN <- "Chile_comun.R"   # ruta al archivo común (Chile_config.R en la misma carpeta)
source(file.path(dirname(PATH_COMUN), "Chile_config.R"), encoding = "UTF-8")
source(file.path(dirname(PATH_COMUN), "Chile_unidades.R"), encoding = "UTF-8")
t0 <- Sys.time()

DIR_HU <- file.path(DIR_OUT, "huerfanos")
dir.create(DIR_HU, recursive = TRUE, showWarnings = FALSE)
f_in <- c(rcri = paste0(TAG_RA, "_EA_1km.tif"), status = paste0(TAG_RA, "_status_EA_1km.tif"),
          red = sprintf("redundancia_AP_%dkm_P%d_EA_1km.tif", UNIT_KM, P))
f_in <- setNames(file.path(DIR_OUT, f_in), names(f_in))
if (!all(file.exists(f_in))) stop("Faltan salidas de APCRS_BCRS_Chile.R en ", DIR_OUT, ".")


# 1 Capas en la rejilla de 1 km ----

st  <- rast(f_in["status"]); rc <- rast(f_in["rcri"]); red <- rast(f_in["red"])
stopifnot(compareGeom(st, rc), compareGeom(st, red))
cell <- which(!is.na(values(st, mat = FALSE)))
v_st  <- values(st,  mat = FALSE)[cell]
v_rc  <- values(rc,  mat = FALSE)[cell]
v_red <- values(red, mat = FALSE)[cell]
stopifnot(identical(v_red == 0, v_st == 1))   # huérfano = redundancia 0

# Clima por vecino más próximo, igual que el resto de capas.
clim <- mask(project(rast(CLIM_VIF), st, method = "near"), st)
X    <- values(clim, mat = TRUE)[cell, , drop = FALSE]
VARS <- names(clim)

# Región, macrozona y AP de cada celda por su centro.
u_rm <- lee_regiones(); reg <- u_rm$reg; mz <- u_rm$mz
aps <- read_sf(PATH_AP_ESC)     |> st_transform(CRS_EA)
f_ap_csv <- file.path(DIR_OUT, paste0(TAG_RA, "_AP.csv"))
if (!setequal(aps[[COL_ID]], read.csv(f_ap_csv)[[COL_ID]]))
  stop(basename(PATH_AP_ESC), " no tiene las mismas AP que la ejecución principal.")
a_reg <- asigna_unidades(reg, st, cell); a_mz <- asigna_unidades(mz, st, cell)
r_reg <- a_reg$r
v_reg <- r_reg[cell]
v_mz  <- a_mz$r[cell]
message("Escenario ", if (BUFFER_KM > 0) paste0("con anillo de ", BUFFER_KM, " km") else "red actual",
        " | superficie protegida (celdas de 1 km) se muestra abajo")
message("Celdas asignadas a la región más próxima ", a_reg$n_fuera, " | a la macrozona más próxima ", a_mz$n_fuera)
v_ap  <- values(rasterize(vect(aps), st, field = 1, background = 0, touches = FALSE), mat = FALSE)[cell] == 1
ex_ap <- terra::extract(st, vect(aps), cells = TRUE, ID = TRUE)
ap_reg <- unique(data.frame(ap = ex_ap$ID, reg = r_reg[ex_ap$cell]))
ap_reg <- ap_reg[!is.na(ap_reg$reg), ]
message(sprintf("Celdas %d | protegidas %d | huérfanas %d", length(cell), sum(v_ap), sum(v_st == 1)))
stopifnot(!anyNA(v_reg), !anyNA(v_mz))


# 2 Tablas por región y macrozona ----

resume <- function(g, nombres, n_ap) {
  n_huer_tot <- sum(v_st == 1)
  do.call(rbind, lapply(seq_along(nombres), function(i) {
    s  <- which(g == i)
    ok <- s[v_st[s] == 0]
    data.frame(unidad = nombres[i], km2 = length(s),
               km2_protegidos = sum(v_ap[s]), pct_protegido = 100 * mean(v_ap[s]), n_AP = n_ap[i],
               km2_huerfanos = sum(v_st[s] == 1), pct_huerfanos = 100 * mean(v_st[s] == 1),
               pct_del_total_huerfanos = 100 * sum(v_st[s] == 1) / n_huer_tot,
               huerfanos_dentro_AP = sum(v_st[s] == 1 & v_ap[s]),
               pct_una_AP = 100 * mean(v_red[s] == 1), pct_cinco_o_mas_AP = 100 * mean(v_red[s] >= 5),
               redundancia_mediana = median(v_red[s]),
               RCRI_mediana = if (length(ok)) median(v_rc[ok]) else NA_real_,
               pct_RCRI_positivo = if (length(ok)) 100 * mean(v_rc[ok] > 0) else NA_real_,
               pct_muy_infrarrep = 100 * mean(v_st[s] == 0 & v_rc[s] < -1))
  }))
}
n_ap_reg <- tabulate(ap_reg$reg, nbins = nrow(reg))
reg_mz   <- vapply(seq_len(nrow(reg)), function(i) {
  x <- v_mz[v_reg == i & !is.na(v_reg)]; if (!length(x)) NA_integer_ else as.integer(names(which.max(table(x))))
}, integer(1))
n_ap_mz  <- vapply(seq_len(nrow(mz)), function(m) length(unique(ap_reg$ap[reg_mz[ap_reg$reg] %in% m])), integer(1))

t_reg <- resume(v_reg, reg[[COL_REGION]], n_ap_reg)
t_reg$macrozona <- mz[[COL_MACROZONA]][reg_mz]
t_mz  <- resume(v_mz, mz[[COL_MACROZONA]], n_ap_mz)
nac   <- resume(rep(1L, length(cell)), "Chile", nrow(aps))
t_reg <- t_reg[order(match(t_reg$macrozona, c("Norte", "Centro", "Centro Sur", "Sur", "Austral")), -t_reg$km2), ]
t_mz  <- rbind(t_mz[order(match(t_mz$unidad, c("Norte", "Centro", "Centro Sur", "Sur", "Austral"))), ], nac)
names(t_reg)[1] <- "region"; names(t_mz)[1] <- "macrozona"
write.csv(t_reg, file.path(DIR_HU, sprintf("resumen_regiones_P%d_%dkm%s.csv", P, UNIT_KM, TAG_ESC)), row.names = FALSE)
write.csv(t_mz,  file.path(DIR_HU, sprintf("resumen_macrozonas_P%d_%dkm%s.csv", P, UNIT_KM, TAG_ESC)), row.names = FALSE)
op <- options(width = 200); print(t_mz, digits = 3, row.names = FALSE); options(op)


# 3 Caracterización climática de los huérfanos ----
# Percentiles 10, 50 y 90 de cada variable en todo Chile, en los climas con
# análogo en la red y en los huérfanos, y de los huérfanos por macrozona. Las
# unidades son las del archivo de clima.

pc <- function(x) quantile(x, c(0.1, 0.5, 0.9), na.rm = TRUE, names = FALSE)
grupos <- list(Chile = rep(TRUE, length(cell)), con_analogo = v_st != 1, huerfanos = v_st == 1)
for (m in seq_len(nrow(mz))) grupos[[paste0("huerfanos_", mz[[COL_MACROZONA]][m])]] <- v_st == 1 & v_mz %in% m
t_cl <- do.call(rbind, lapply(names(grupos), function(g) {
  s <- grupos[[g]]
  if (!any(s, na.rm = TRUE)) return(NULL)
  do.call(rbind, lapply(VARS, function(v) {
    q <- pc(X[s, v]); data.frame(grupo = g, celdas = sum(s, na.rm = TRUE), variable = v, p10 = q[1], p50 = q[2], p90 = q[3])
  }))
}))
write.csv(t_cl, file.path(DIR_HU, sprintf("clima_huerfanos_P%d_%dkm%s.csv", P, UNIT_KM, TAG_ESC)), row.names = FALSE)
print(reshape(t_cl[t_cl$grupo %in% c("Chile", "con_analogo", "huerfanos"), c("grupo", "variable", "p50")],
              idvar = "variable", timevar = "grupo", direction = "wide"), digits = 4, row.names = FALSE)


# 4 Capa de clases de redundancia ----

clase <- cut(v_red, c(-Inf, 0, 1, 4, Inf), labels = FALSE)   # 1 huérfano, 2 una AP, 3 de 2 a 4, 4 cinco o más
r_cl  <- rast(st); vv <- rep(NA_real_, ncell(r_cl)); vv[cell] <- clase; values(r_cl) <- vv
levels(r_cl) <- data.frame(id = 1:4, clase = c("Huérfano", "Una AP", "2 a 4 AP", "5 o más AP"))
col_cl <- c("#B2182B", "#EF8A62", "#92C5DE", "#2166AC")
coltab(r_cl) <- data.frame(value = 1:4, col = col_cl)
writeRaster(r_cl, file.path(DIR_HU, sprintf("clases_redundancia_P%d_%dkm%s_EA_1km.tif", P, UNIT_KM, TAG_ESC)),
            datatype = "INT1U", overwrite = TRUE, gdal = "COMPRESS=LZW")


# 5 Tipos de clima de Köppen ----
# Cada celda toma el tipo de Köppen de su centro (las de borde, el del polígono
# más próximo). Los códigos se simplifican a la letra principal y dos más
# (BWk(s) pasa a BWk).

v_kp <- NULL
kp <- lee_koppen()
if (!is.null(kp)) {
  v_kp  <- asigna_unidades(kp, st, cell)$r[cell]
  tipos <- kp$tipo
  desc  <- kp$descripcion

  t_kp <- resume(v_kp, tipos, rep(NA_integer_, length(tipos)))
  t_kp$n_AP <- NULL
  t_kp <- cbind(tipo = t_kp$unidad, descripcion = desc, t_kp[, -1])
  t_kp <- t_kp[order(-t_kp$km2_huerfanos, -t_kp$km2), ]
  write.csv(t_kp, file.path(DIR_HU, sprintf("resumen_koppen_P%d_%dkm%s.csv", P, UNIT_KM, TAG_ESC)), row.names = FALSE)
  op <- options(width = 200)
  print(t_kp[, c("tipo", "descripcion", "km2", "pct_protegido", "km2_huerfanos", "pct_huerfanos",
                 "pct_del_total_huerfanos", "pct_muy_infrarrep", "RCRI_mediana")], digits = 3, row.names = FALSE)
  options(op)

  # Colores. Tipos con al menos el 2 % de Chile o el 1 % de los huérfanos; el
  # resto, en gris como «Otros».
  imp <- t_kp$tipo[t_kp$km2 >= 0.02 * length(cell) | t_kp$pct_del_total_huerfanos >= 1]
  imp <- tipos[tipos %in% imp]
  imp <- imp[order(substr(imp, 1, 1), imp)]
  # Colores por grupo de Köppen, para que el gráfico se lea por familias de clima.
  grupos_col <- list(A = c("#B2182B", "#D6604D"),
                     B = c("#8C510A", "#E08214", "#DFC27D", "#BF812D", "#FDB863", "#7F3B08"),
                     C = c("#1B7837", "#5AAE61", "#A6DBA0", "#01665E", "#35978F", "#80CDC1", "#276419"),
                     D = c("#762A83", "#9970AB", "#C2A5CF", "#40004B"),
                     E = c("#2166AC", "#92C5DE", "#053061"))
  pal <- unlist(lapply(split(imp, substr(imp, 1, 1)), function(g) {
    cc <- grupos_col[[substr(g[1], 1, 1)]]
    if (is.null(cc)) cc <- "grey40"
    setNames(rep(cc, length.out = length(g)), g)
  }), use.names = FALSE)
  names(pal) <- unlist(split(imp, substr(imp, 1, 1)), use.names = FALSE)
  pal <- pal[imp]
  col_kp <- ifelse(tipos[v_kp] %in% imp, pal[tipos[v_kp]], "grey70")
  leyenda <- c(paste0(imp, "  ", ifelse(imp %in% names(NOMBRES_KOPPEN), NOMBRES_KOPPEN[imp], "")), "Otros")
  col_ley <- c(pal, "grey70")
}


# 6 Figuras ----

png(file.path(DIR_HU, sprintf("mapa_clases_redundancia_P%d_%dkm%s.png", P, UNIT_KM, TAG_ESC)),
    width = 1600, height = 3200, res = 300)
plot(r_cl, axes = FALSE, box = FALSE, mar = c(1, 1, 1, 6),
     plg = list(title = "Climas respaldados por", cex = 0.7))
plot(st_geometry(reg), add = TRUE, border = "grey30", lwd = 0.4)
plot(st_geometry(mz), add = TRUE, border = "black", lwd = 0.9)
dev.off()

ETIQ <- c(bio1 = "Temperatura media anual (°C)", bio2 = "Oscilación térmica diaria media (°C)",
          bio4 = "Estacionalidad de la temperatura (DE × 100)", bio12 = "Precipitación anual (mm)",
          bio15 = "Estacionalidad de la precipitación (CV, %)")
pares <- Filter(function(p) all(p %in% VARS), list(c("bio1", "bio12"), c("bio4", "bio15"), c("bio1", "bio2")))
set.seed(1)
todo  <- sample(seq_along(cell), min(200000, length(cell)))
fondo <- sample(which(v_st != 1), min(150000, sum(v_st != 1)))
hu    <- sample(which(v_st == 1))
panel <- function(p, i_fondo, col_fondo, i_frente = NULL, col_frente = NULL) {
  plot(X[i_fondo, p[1]], X[i_fondo, p[2]], pch = ".", col = col_fondo, xlab = ETIQ[p[1]], ylab = ETIQ[p[2]],
       xlim = range(X[, p[1]], na.rm = TRUE), ylim = range(X[, p[2]], na.rm = TRUE))
  if (!is.null(i_frente)) points(X[i_frente, p[1]], X[i_frente, p[2]], pch = ".", col = col_frente)
}
f_esp <- file.path(DIR_HU, sprintf("espacio_climatico_huerfanos_P%d_%dkm%s.png", P, UNIT_KM, TAG_ESC))
if (is.null(v_kp)) {
  png(f_esp, width = 1300 * length(pares), height = 1300, res = 300)
  par(mfrow = c(1, length(pares)), mar = c(4.2, 4.2, 1, 1), mgp = c(2.4, 0.7, 0), cex.lab = 0.8)
  for (p in pares) panel(p, fondo, adjustcolor("grey55", 0.4), hu, adjustcolor(col_cl[1], 0.5))
  dev.off()
} else {
  png(f_esp, width = 1300 * length(pares), height = 2900, res = 300)
  layout(rbind(matrix(seq_len(2 * length(pares)), 2, byrow = TRUE), 2 * length(pares) + 1),
         heights = c(1, 1, 0.3))
  par(mar = c(4.2, 4.2, 2, 1), mgp = c(2.4, 0.7, 0), cex.lab = 0.8)
  for (k in seq_along(pares)) {
    panel(pares[[k]], todo, adjustcolor(col_kp[todo], 0.5))
    if (k == 1) mtext("Todo Chile por tipo de clima de Köppen", side = 3, adj = 0, line = 0.5, cex = 0.75, font = 2)
  }
  for (k in seq_along(pares)) {
    panel(pares[[k]], fondo, adjustcolor("grey80", 0.4), hu, adjustcolor(col_kp[hu], 0.7))
    if (k == 1) mtext("Climas sin protección (huérfanos) por tipo de clima de Köppen", side = 3, adj = 0,
                      line = 0.5, cex = 0.75, font = 2)
  }
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", legend = leyenda, fill = col_ley, border = NA, bty = "n", ncol = 3, cex = 1.1)
  dev.off()

  # Barras por tipo de Köppen, ordenadas por la fracción sin protección suficiente.
  tb <- t_kp[t_kp$km2 >= 0.005 * length(cell), ]
  tb$pct_resto <- 100 - tb$pct_huerfanos - tb$pct_muy_infrarrep
  tb <- tb[order(tb$pct_huerfanos + tb$pct_muy_infrarrep), ]
  m  <- t(as.matrix(tb[, c("pct_huerfanos", "pct_muy_infrarrep", "pct_resto")]))
  png(file.path(DIR_HU, sprintf("koppen_huerfanos_P%d_%dkm%s.png", P, UNIT_KM, TAG_ESC)),
      width = 2600, height = 500 + 120 * nrow(tb), res = 300)
  layout(matrix(1:2, 2), heights = c(120 * nrow(tb) + 300, 200))
  par(mar = c(4, 17, 1, 1), mgp = c(2.2, 0.6, 0))
  barplot(m, horiz = TRUE, names.arg = paste0(tb$descripcion, " (", tb$tipo, ")"), las = 1, cex.names = 0.7,
          col = c(col_cl[1], col_cl[2], "grey85"), border = NA, xlab = "% de la superficie del tipo de clima",
          xlim = c(0, 100))
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", c("Sin análogo en la red", "Menos de 1/10 de la protección esperada", "Resto"),
         fill = c(col_cl[1], col_cl[2], "grey85"), border = NA, bty = "n", cex = 0.7, horiz = TRUE)
  dev.off()
}

message("Terminado | ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
