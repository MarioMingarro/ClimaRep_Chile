# Regiones, macrozonas y tipos de Köppen en la rejilla de 1 km. Lo cargan los
# scripts que resumen por unidades (Huerfanos_regiones_Chile.R y
# Comparar_escenarios_Chile.R), para que usen exactamente la misma asignación.
# Requiere Chile_config.R cargado.

NOMBRES_KOPPEN <- c(
  Af = "Tropical húmedo", Am = "Monzónico", Aw = "Sabana",
  BWh = "Desierto cálido", BWk = "Desierto frío", BWn = "Desierto costero nuboso",
  BSh = "Estepa cálida", BSk = "Estepa fría", BSn = "Estepa costera nubosa",
  Csa = "Mediterráneo de verano cálido", Csb = "Mediterráneo de verano templado",
  Csc = "Mediterráneo de verano frío", Cwa = "Subtropical de invierno seco",
  Cwb = "Templado de altura de invierno seco", Cwc = "Frío de altura de invierno seco",
  Cfa = "Subtropical húmedo", Cfb = "Oceánico templado", Cfc = "Oceánico subpolar",
  Dsa = "Continental de verano seco y cálido", Dsb = "Continental de verano seco y templado",
  Dsc = "Continental de verano seco y frío", Dwb = "Continental de invierno seco",
  Dwc = "Subártico de invierno seco", Dfb = "Continental húmedo", Dfc = "Subártico",
  ET = "Tundra", EF = "Hielo perpetuo")

# Unidad de cada celda por su centro. Las celdas de costa o frontera que caen
# fuera de la capa se asignan al polígono más próximo, para que las tablas
# sumen el total nacional. Devuelve el vector para toda la rejilla.
asigna_unidades <- function(pol, st, cell) {
  r <- values(rasterize(vect(pol), st, field = seq_len(nrow(pol)), touches = FALSE), mat = FALSE)
  fuera <- intersect(which(is.na(r)), cell)
  if (length(fuera)) {
    xy <- st_as_sf(as.data.frame(xyFromCell(st, fuera)), coords = c("x", "y"), crs = CRS_EA)
    r[fuera] <- st_nearest_feature(xy, pol)
  }
  list(r = r, n_fuera = length(fuera))
}

lee_regiones <- function() list(
  reg = read_sf(PATH_REGIONES)   |> st_make_valid() |> st_transform(CRS_EA),
  mz  = read_sf(PATH_MACROZONAS) |> st_make_valid() |> st_transform(CRS_EA))

# Tipos de Köppen simplificados a tres letras (BWk(s) pasa a BWk), un polígono
# por tipo. NULL si no está la capa.
lee_koppen <- function() {
  if (!file.exists(PATH_KOPPEN)) return(NULL)
  kp <- read_sf(PATH_KOPPEN) |> st_make_valid() |> st_transform(CRS_EA)
  col <- COL_KOPPEN
  col_txt <- names(kp)[vapply(st_drop_geometry(kp), is.character, logical(1))]
  if (is.null(col)) {
    parecido <- vapply(col_txt, function(c) mean(grepl("^[A-E][A-Za-z]{0,2}", kp[[c]])), numeric(1))
    col <- if (length(parecido) && max(parecido) > 0.5) names(which.max(parecido)) else col_txt[1]
  }
  cod  <- kp[[col]]
  base <- ifelse(grepl("^[A-E][A-Za-z]{0,2}", cod), regmatches(cod, regexpr("^[A-E][A-Za-z]{0,2}", cod)), cod)
  message("Köppen | campo ", col, " | tipos: ", paste(sort(unique(base)), collapse = ", "))
  kp$tipo <- base
  kp <- aggregate(kp[, "tipo"], by = list(t = kp$tipo), FUN = function(x) x[1])[, "tipo"]
  kp$descripcion <- ifelse(kp$tipo %in% names(NOMBRES_KOPPEN), NOMBRES_KOPPEN[kp$tipo], kp$tipo)
  kp
}
