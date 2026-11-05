comprobar <- function(condicion, mensaje) {
  if (!isTRUE(condicion)) stop(mensaje, call. = FALSE)
}

validar_manifiesto <- function(x) {
  columnas <- c("archivo", "grupo", "concentracion_g_l",
                "densidad_g_cm3", "repeticion")
  comprobar(all(columnas %in% names(x)), "Manifiesto: faltan columnas.")
  comprobar(nrow(x) == 25L && !anyNA(x), "Se requieren 25 ensayos completos.")
  comprobar(!anyDuplicated(x$archivo), "Hay archivos duplicados.")
  comprobar(all(grepl("^imu_[0-9]{8}_[0-9]{6}\\.csv$", x$archivo)),
            "Nombre de archivo no permitido.")
  comprobar(identical(sort(unique(x$grupo)), paste0("G", 1:5)),
            "Se requieren los grupos G1 a G5.")
  for (g in paste0("G", 1:5)) {
    d <- x[x$grupo == g, ]
    comprobar(nrow(d) == 5L && identical(sort(d$repeticion), 1:5),
              paste("Repeticiones incorrectas:", g))
    comprobar(length(unique(d$densidad_g_cm3)) == 1L &&
                length(unique(d$concentracion_g_l)) == 1L,
              paste("Condiciones inconsistentes:", g))
  }
  comprobar(all(is.finite(x$densidad_g_cm3)) && all(x$densidad_g_cm3 > 0) &&
              length(unique(x$densidad_g_cm3)) == 5L,
            "Densidades no validas.")
}

preparar_registro <- function(d, nombre) {
  requeridas <- c("t_ms", "ax", "ay", "az", "gx", "gy", "gz", "temp", "ang_deg")
  comprobar(all(requeridas %in% names(d)), paste(nombre, ": faltan columnas."))
  comprobar(nrow(d) > 2L, paste(nombre, ": registro insuficiente."))
  comprobar(all(vapply(d[requeridas], is.numeric, logical(1))),
            paste(nombre, ": columnas no numericas."))
  comprobar(all(vapply(d[requeridas], function(x) all(is.finite(x)), logical(1))),
            paste(nombre, ": hay valores ausentes o no finitos; revisar el origen."))
  comprobar(all(d$t_ms >= 0), paste(nombre, ": tiempos negativos."))
  dt <- diff(d$t_ms)
  segmento <- cumsum(c(TRUE, dt <= 0))
  ultimo <- d[segmento == max(segmento), ]
  ultimo$t_s <- (ultimo$t_ms - ultimo$t_ms[1]) / 1000
  comprobar(max(ultimo$t_s) >= 30, paste(nombre, ": ultimo segmento menor de 30 s."))
  list(
    datos = ultimo,
    calidad = data.frame(
      archivo = nombre, muestras_originales = nrow(d),
      reinicios_o_duplicados = sum(dt <= 0),
      segmentos = max(segmento), muestras_segmento_final = nrow(ultimo),
      muestras_previas_no_usadas = nrow(d) - nrow(ultimo),
      duracion_segmento_final_s = max(ultimo$t_s),
      mediana_dt_ms = median(diff(ultimo$t_ms)),
      huecos_mayores_50ms = sum(dt > 50),
      huecos_segmento_final_mayores_50ms = sum(diff(ultimo$t_ms) > 50)
    )
  )
}

resumir_ensayo <- function(d, ventanas = c(5, 10, 15, 20)) {
  comprobar(all(diff(d$t_s) > 0), "El tiempo del segmento no es creciente.")
  fin <- max(d$t_s)
  dt_nominal <- median(diff(d$t_s))
  medias <- vapply(ventanas, function(w) {
    comprobar(fin >= w, "Ventana mayor que la duracion del registro.")
    z <- d[d$t_s > fin - w, ]
    comprobar(nrow(z) >= 10L && nrow(z) * dt_nominal / w >= 0.9 &&
                max(diff(z$t_s)) <= 1 &&
                min(z$t_s) - (fin - w) <= 1,
              paste("Cobertura insuficiente en ventana terminal de", w, "s."))
    mean(z$ang_deg)
  }, numeric(1))
  hueco_maximo <- max(vapply(ventanas, function(w) {
    z <- d[d$t_s > fin - w, ]
    max(c(diff(z$t_s), min(z$t_s) - (fin - w)))
  }, numeric(1)))
  cobertura <- hueco_maximo <= 0.1 + 1e-8
  terminal <- d[d$t_s > fin - 10, ]
  pendiente <- unname(coef(lm(ang_deg ~ t_s, data = terminal))[2])
  diferencia <- max(abs(medias - medias[ventanas == 10]))
  estable <- abs(pendiente * 10) <= 0.1 && diferencia <= 0.1 && cobertura
  salida <- data.frame(
    angulo_terminal_deg = medias[ventanas == 10],
    sd_intraensayo_deg = sd(terminal$ang_deg),
    muestras_terminales = nrow(terminal),
    pendiente_terminal_deg_s = pendiente,
    deriva_lineal_10s_deg = pendiente * 10,
    diferencia_max_ventanas_deg = diferencia,
    max_hueco_ventanas_s = hueco_maximo,
    cobertura_orientativa = cobertura,
    estabilidad_orientativa = estable,
    temp_sensor_terminal_media = mean(terminal$temp)
  )
  for (i in seq_along(ventanas)) salida[[paste0("media_", ventanas[i], "s_deg")]] <- medias[i]
  salida
}

resumir_grupos <- function(d, variable = "angulo_terminal_deg") {
  do.call(rbind, lapply(split(d, d$grupo), function(z) {
    x <- z[[variable]]
    n <- length(x)
    sem <- sd(x) / sqrt(n)
    margen <- qt(0.975, df = n - 1) * sem
    data.frame(
      grupo = z$grupo[1], densidad_g_cm3 = z$densidad_g_cm3[1], n = n,
      media_deg = mean(x), sd_entre_ensayos_deg = sd(x), sem_deg = sem,
      ic95_inferior_deg = mean(x) - margen,
      ic95_superior_deg = mean(x) + margen
    )
  }))
}

comparar_pares <- function(d) {
  pares <- combn(sort(unique(d$grupo)), 2, simplify = FALSE)
  resultado <- do.call(rbind, lapply(pares, function(par) {
    x <- d$angulo_terminal_deg[d$grupo == par[1]]
    y <- d$angulo_terminal_deg[d$grupo == par[2]]
    prueba <- t.test(x, y, var.equal = FALSE)
    data.frame(
      grupo_a = par[1], grupo_b = par[2],
      diferencia_a_menos_b_deg = mean(x) - mean(y),
      ic95_nominal_inferior_deg = prueba$conf.int[1],
      ic95_nominal_superior_deg = prueba$conf.int[2],
      t = unname(prueba$statistic), gl = unname(prueba$parameter),
      p_sin_ajustar = prueba$p.value
    )
  }))
  resultado$p_holm <- p.adjust(resultado$p_sin_ajustar, method = "holm")
  resultado
}

validar_calibracion <- function(d, modo) {
  comprobar(modo %in% c("ensayo", "densidad"), "Modo de validacion desconocido.")
  do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
    excluir <- if (modo == "ensayo") seq_len(nrow(d)) == i else d$grupo == d$grupo[i]
    entrenamiento <- d[!excluir, ]
    modelo <- lm(densidad_g_cm3 ~ angulo_terminal_deg, data = entrenamiento)
    pred <- unname(predict(modelo, newdata = d[i, ]))
    data.frame(
      archivo = d$archivo[i], grupo = d$grupo[i], modo = modo,
      n_entrenamiento = nrow(entrenamiento),
      niveles_entrenamiento = length(unique(entrenamiento$grupo)),
      densidad_referencia_g_cm3 = d$densidad_g_cm3[i],
      densidad_predicha_g_cm3 = pred,
      error_g_cm3 = pred - d$densidad_g_cm3[i],
      referencia_fuera_rango_entrenamiento =
        d$densidad_g_cm3[i] < min(entrenamiento$densidad_g_cm3) ||
        d$densidad_g_cm3[i] > max(entrenamiento$densidad_g_cm3)
    )
  }))
}

pdf_texto <- function(lineas, ruta) {
  envueltas <- unlist(lapply(lineas, function(x) {
    if (!nzchar(x)) "" else strwrap(x, width = 94)
  }), use.names = FALSE)
  pdf(ruta, width = 8.27, height = 11.69, family = "Courier")
  on.exit(dev.off())
  paginas <- split(envueltas, ceiling(seq_along(envueltas) / 58))
  for (i in seq_along(paginas)) {
    par(mar = c(2, 2, 2, 2))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    text(0, seq(0.98, 0.03, length.out = 58)[seq_along(paginas[[i]])],
         paginas[[i]], adj = c(0, 1), cex = 0.7)
    mtext(paste("Borrador - Modulo III |", i, "/", length(paginas)), side = 1, cex = 0.7)
  }
}
