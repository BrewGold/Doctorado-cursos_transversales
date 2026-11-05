#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = FALSE)
archivo <- sub("^--file=", "", args[grepl("^--file=", args)])
raiz <- dirname(dirname(normalizePath(archivo)))
source(file.path(raiz, "scripts", "funciones.R"))

esperar_error <- function(expr) {
  resultado <- tryCatch({ force(expr); FALSE }, error = function(e) TRUE)
  stopifnot(resultado)
}
cerca <- function(x, y, tolerancia = 1e-9) {
  stopifnot(isTRUE(all.equal(x, y, tolerance = tolerancia, check.attributes = FALSE)))
}
registro <- function(t, angulo = rep(60, length(t))) {
  data.frame(t_ms = t, ax = 0, ay = 0, az = 9.81, gx = 0, gy = 0,
             gz = 0, temp = 30, ang_deg = angulo)
}

manifiesto <- read.csv(file.path(raiz, "datos", "ensayos.csv"))
validar_manifiesto(manifiesto)
malo <- manifiesto
malo$archivo[2] <- malo$archivo[1]
esperar_error(validar_manifiesto(malo))
malo <- manifiesto
malo$repeticion[1] <- 2
esperar_error(validar_manifiesto(malo))
malo <- manifiesto
malo$densidad_g_cm3[1] <- NA
esperar_error(validar_manifiesto(malo))

constante <- registro(seq(0, 40000, by = 10))
preparado <- preparar_registro(constante, "constante")
stopifnot(preparado$calidad$segmentos == 1, preparado$calidad$mediana_dt_ms == 10)
s <- resumir_ensayo(preparado$datos)
cerca(s$angulo_terminal_deg, 60)
cerca(s$sd_intraensayo_deg, 0)
stopifnot(s$muestras_terminales == 1000, s$estabilidad_orientativa)

reinicio <- preparar_registro(rbind(registro(seq(0, 1000, 10), rep(999, 101)),
                                   constante), "reinicio")
stopifnot(reinicio$calidad$reinicios_o_duplicados == 1,
          reinicio$calidad$muestras_previas_no_usadas == 101)
cerca(resumir_ensayo(reinicio$datos)$angulo_terminal_deg, 60)
duplicado <- constante[c(1:100, 100, 101:nrow(constante)), ]
stopifnot(preparar_registro(duplicado, "duplicado")$calidad$segmentos == 2)

malo <- constante
malo$ang_deg[1] <- NA
esperar_error(preparar_registro(malo, "NA"))
malo$ang_deg[1] <- Inf
esperar_error(preparar_registro(malo, "Inf"))
esperar_error(preparar_registro(constante[, -1], "sin_tiempo"))
esperar_error(preparar_registro(constante[1:100, ], "corto"))

tendencia <- preparar_registro(registro(seq(0, 40000, 10),
                                       60 + seq(0, 40000, 10) / 1000), "deriva")
stopifnot(!resumir_ensayo(tendencia$datos)$estabilidad_orientativa)
cerca(resumir_ensayo(tendencia$datos)$deriva_lineal_10s_deg, 10)
hueco <- constante[!(constante$t_ms > 35000 & constante$t_ms < 35140), ]
h <- resumir_ensayo(preparar_registro(hueco, "hueco")$datos)
stopifnot(!h$cobertura_orientativa, !h$estabilidad_orientativa)
cerca(h$max_hueco_ventanas_s, 0.14)
cerca(h$angulo_terminal_deg, 60)
hueco_grande <- constante[!(constante$t_ms > 35000 & constante$t_ms < 37000), ]
esperar_error(resumir_ensayo(preparar_registro(hueco_grande, "hueco_grande")$datos))
escaso <- constante[!(constante$t_ms > 30000 & seq_len(nrow(constante)) %% 2 == 0), ]
esperar_error(resumir_ensayo(preparar_registro(escaso, "escaso")$datos))

sintetico <- manifiesto
sintetico$angulo_terminal_deg <- rep(1:5, each = 5) + rep(c(-0.2, -0.1, 0, 0.1, 0.2), 5)
r <- resumir_grupos(sintetico)
cerca(r$media_deg, 1:5)
cerca(r$sd_entre_ensayos_deg, rep(sd(c(-0.2, -0.1, 0, 0.1, 0.2)), 5))
cerca(r$ic95_superior_deg - r$media_deg, qt(0.975, 4) * r$sem_deg)
pares <- comparar_pares(sintetico)
stopifnot(nrow(pares) == 10, all(pares$p_holm >= pares$p_sin_ajustar))
cerca(pares$p_holm, p.adjust(pares$p_sin_ajustar, "holm"))

for (modo in c("ensayo", "densidad")) {
  cv <- validar_calibracion(sintetico, modo)
  stopifnot(nrow(cv) == 25,
            all(cv$n_entrenamiento == if (modo == "ensayo") 24 else 20),
            all(cv$niveles_entrenamiento == if (modo == "ensayo") 5 else 4))
  for (i in seq_len(nrow(sintetico))) {
    quitar <- if (modo == "ensayo") seq_len(nrow(sintetico)) == i else
      sintetico$grupo == sintetico$grupo[i]
    ajuste <- lm(densidad_g_cm3 ~ angulo_terminal_deg, data = sintetico[!quitar, ])
    cerca(cv$densidad_predicha_g_cm3[i], unname(predict(ajuste, sintetico[i, ])))
    cambiado <- sintetico
    cambiado$densidad_g_cm3[quitar] <- 999
    sin_fuga <- validar_calibracion(cambiado, modo)
    cerca(sin_fuga$densidad_predicha_g_cm3[i], cv$densidad_predicha_g_cm3[i])
  }
}

salidas <- file.path(raiz, "resultados")
leer <- function(nombre) read.csv(file.path(salidas, nombre), stringsAsFactors = FALSE)
e <- leer("ensayos.csv")
q <- leer("calidad.csv")
g <- leer("resumen_grupos.csv")
stopifnot(nrow(e) == 25, identical(e$archivo, manifiesto$archivo),
          nrow(q) == 25, nrow(g) == 5, all(g$n == 5),
          sum(q$reinicios_o_duplicados) == 2,
          sum(!e$cobertura_orientativa) == 3,
          all(is.finite(e$angulo_terminal_deg)))
for (i in 1:5) {
  z <- e$angulo_terminal_deg[e$grupo == paste0("G", i)]
  cerca(g$media_deg[i], mean(z))
  cerca(g$sd_entre_ensayos_deg[i], sd(z))
}
cerca(e$angulo_terminal_deg, e$media_10s_deg)
cerca(e$angulo_terminal_deg[1], 63.1152679486468)
global <- leer("contrastes_globales.csv")
cerca(global$p[1], oneway.test(angulo_terminal_deg ~ factor(grupo),
                               data = e, var.equal = FALSE)$p.value)
stopifnot(nrow(leer("comparaciones_pares.csv")) == 10,
          nrow(leer("normalidad_por_grupo.csv")) == 5,
          nrow(leer("sensibilidad_ventanas.csv")) == 4,
          nrow(leer("calibracion_validacion.csv")) == 50)
cv <- leer("calibracion_validacion.csv")
metricas <- leer("calibracion_metricas.csv")
for (i in seq_len(nrow(metricas))) {
  errores <- cv$error_g_cm3[cv$modo == metricas$modo[i]]
  cerca(metricas$rmse_g_cm3[i], sqrt(mean(errores^2)))
}
memoria <- readLines(file.path(salidas, "memoria_borrador.md"))
stopifnot(any(grepl("Anexo A. Codigo R completo", memoria, fixed = TRUE)),
          any(grepl("### scripts/analizar.R", memoria, fixed = TRUE)),
          any(grepl("### scripts/funciones.R", memoria, fixed = TRUE)))
for (nombre in c("figuras.pdf", "memoria_borrador.pdf", "propuesta_borrador.pdf")) {
  ruta <- file.path(salidas, nombre)
  stopifnot(file.exists(ruta), file.info(ruta)$size > 1000)
  con <- file(ruta, "rb")
  cabecera <- rawToChar(readBin(con, "raw", 5))
  close(con)
  stopifnot(cabecera == "%PDF-")
}
message("Pruebas correctas: datos sinteticos, errores, ausencia de fuga y resultados completos.")
