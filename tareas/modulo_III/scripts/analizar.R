#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = FALSE)
archivo_script <- sub("^--file=", "", args[grepl("^--file=", args)])
if (length(archivo_script) != 1L) stop("Ejecutar con Rscript scripts/analizar.R")
raiz <- dirname(dirname(normalizePath(archivo_script)))
source(file.path(raiz, "scripts", "funciones.R"))
options(warn = 1)

revision <- "bcc97816212b48ea0aed2849c85d9b0065c223ed"
url <- paste0("https://raw.githubusercontent.com/BrewGold/Doctorado-densimetro/",
              revision, "/references/tfm/datos_udp.zip")
md5_esperado <- "2aeb5372c77fe3b04d520eade81e3b18"
cache <- file.path(raiz, "datos", "cache")
salidas <- file.path(raiz, "resultados")
dir.create(cache, recursive = TRUE, showWarnings = FALSE)
dir.create(salidas, recursive = TRUE, showWarnings = FALSE)
zip <- file.path(cache, "datos_udp.zip")
if (!file.exists(zip)) {
  message("Descargando archivo experimental de la revision fijada...")
  parcial <- paste0(zip, ".part")
  estado <- download.file(url, parcial, mode = "wb", quiet = TRUE)
  comprobar(estado == 0L, "Fallo de descarga; revisar conexion y archivo .part.")
  comprobar(identical(unname(tools::md5sum(parcial)), md5_esperado),
            "El archivo descargado no coincide con el MD5 esperado.")
  comprobar(file.rename(parcial, zip), "No se pudo guardar la descarga.")
}
comprobar(identical(unname(tools::md5sum(zip)), md5_esperado),
          "La cache no coincide con el MD5 esperado; revisar antes de continuar.")
manifiesto <- read.csv(file.path(raiz, "datos", "ensayos.csv"), stringsAsFactors = FALSE)
validar_manifiesto(manifiesto)
miembros <- unzip(zip, list = TRUE)$Name
esperados <- paste0("datos_udp/", manifiesto$archivo)
comprobar(setequal(miembros[grepl("\\.csv$", miembros)], esperados),
          "Los CSV del archivo no coinciden con el manifiesto.")
unzip(zip, files = esperados, exdir = cache, junkpaths = TRUE)
comprobar(all(file.exists(file.path(cache, manifiesto$archivo))),
          "No se pudieron extraer todos los ensayos.")
escribir <- function(d, nombre) {
  write.csv(d, file.path(salidas, nombre), row.names = FALSE, na = "")
}

registros <- lapply(seq_len(nrow(manifiesto)), function(i) {
  preparar_registro(read.csv(file.path(cache, manifiesto$archivo[i]),
                              stringsAsFactors = FALSE), manifiesto$archivo[i])
})
calidad <- do.call(rbind, lapply(registros, `[[`, "calidad"))
ensayos <- cbind(manifiesto, do.call(rbind, lapply(registros, function(r) {
  resumir_ensayo(r$datos)
})))
escribir(calidad, "calidad.csv")
escribir(ensayos, "ensayos.csv")
resumen <- resumir_grupos(ensayos)
escribir(resumen, "resumen_grupos.csv")

welch <- oneway.test(angulo_terminal_deg ~ factor(grupo), data = ensayos, var.equal = FALSE)
kw <- kruskal.test(angulo_terminal_deg ~ factor(grupo), data = ensayos)
contrastes <- data.frame(
  prueba = c("ANOVA_Welch", "Kruskal_Wallis_sensibilidad"),
  estadistico = c(unname(welch$statistic), unname(kw$statistic)),
  gl1 = c(unname(welch$parameter[1]), unname(kw$parameter)),
  gl2 = c(unname(welch$parameter[2]), NA_real_),
  p = c(welch$p.value, kw$p.value)
)
escribir(contrastes, "contrastes_globales.csv")
pares <- comparar_pares(ensayos)
escribir(pares, "comparaciones_pares.csv")
normalidad <- do.call(rbind, lapply(split(ensayos, ensayos$grupo), function(z) {
  s <- shapiro.test(z$angulo_terminal_deg)
  data.frame(grupo = z$grupo[1], n = nrow(z), W = unname(s$statistic), p = s$p.value)
}))
escribir(normalidad, "normalidad_por_grupo.csv")
ensayos$desviacion_mediana <- abs(ensayos$angulo_terminal_deg -
  ave(ensayos$angulo_terminal_deg, ensayos$grupo, FUN = median))
bf <- anova(lm(desviacion_mediana ~ factor(grupo), data = ensayos))
escribir(data.frame(prueba = "Brown_Forsythe", F = bf$`F value`[1],
                    gl1 = bf$Df[1], gl2 = bf$Df[2], p = bf$`Pr(>F)`[1]),
         "homogeneidad_varianzas.csv")
sensibilidad <- do.call(rbind, lapply(c(5, 10, 15, 20), function(w) {
  x <- ensayos[[paste0("media_", w, "s_deg")]]
  prueba <- oneway.test(x ~ factor(ensayos$grupo), var.equal = FALSE)
  data.frame(ventana_s = w, F_welch = unname(prueba$statistic), p = prueba$p.value,
             cambio_max_vs_10s_deg = max(abs(x - ensayos$angulo_terminal_deg)))
}))
escribir(sensibilidad, "sensibilidad_ventanas.csv")

modelo <- lm(densidad_g_cm3 ~ angulo_terminal_deg, data = ensayos)
coeficientes <- as.data.frame(coef(summary(modelo)))
names(coeficientes) <- c("estimacion", "error_estandar", "t", "p")
coeficientes$termino <- rownames(coeficientes)
coeficientes <- coeficientes[, c("termino", "estimacion", "error_estandar", "t", "p")]
escribir(coeficientes, "calibracion_coeficientes.csv")
cv <- rbind(validar_calibracion(ensayos, "ensayo"),
            validar_calibracion(ensayos, "densidad"))
escribir(cv, "calibracion_validacion.csv")
metricas <- do.call(rbind, lapply(split(cv, cv$modo), function(z) {
  data.frame(modo = z$modo[1], n_predicciones = nrow(z), rmse_g_cm3 = sqrt(mean(z$error_g_cm3^2)),
             mae_g_cm3 = mean(abs(z$error_g_cm3)), sesgo_g_cm3 = mean(z$error_g_cm3))
}))
escribir(metricas, "calibracion_metricas.csv")
cv_grupos <- do.call(rbind, lapply(split(cv[cv$modo == "densidad", ], cv$grupo[cv$modo == "densidad"]),
  function(z) data.frame(
    grupo = z$grupo[1], rmse_g_cm3 = sqrt(mean(z$error_g_cm3^2)),
    sesgo_g_cm3 = mean(z$error_g_cm3),
    referencia_fuera_rango_entrenamiento = all(z$referencia_fuera_rango_entrenamiento)
  )))
escribir(cv_grupos, "calibracion_por_densidad.csv")

graficar <- function() {
  pdf(file.path(salidas, "figuras.pdf"), width = 10, height = 8)
  on.exit(dev.off())
  colores <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00")
  for (g in paste0("G", 1:5)) {
    par(mfrow = c(3, 2), mar = c(3.5, 4, 3, 1))
    for (i in which(manifiesto$grupo == g)) {
      d <- registros[[i]]$datos
      indices <- unique(c(seq(1, nrow(d), length.out = min(1500L, nrow(d))), nrow(d)))
      plot(d$t_s[indices], d$ang_deg[indices], type = "l", col = colores[match(g, paste0("G", 1:5))],
           xlab = "Tiempo del segmento final (s)", ylab = "ang_deg (grados)",
           main = manifiesto$archivo[i], cex.main = 0.8)
      abline(v = max(d$t_s) - 10, lty = 2)
      abline(h = ensayos$angulo_terminal_deg[i], col = "red", lty = 3)
    }
    plot.new()
    text(0.5, 0.6, paste(g, "- ultimo segmento temporal\nLinea vertical: ventana final de 10 s\n",
                        "Linea roja: media terminal\nGrafico reducido; calculos con todas las muestras",
                        "\nNo se interpolan huecos."), cex = 0.9)
  }
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  plot(seq_len(nrow(ensayos)), ensayos$angulo_terminal_deg, pch = 19,
       col = colores[match(ensayos$grupo, paste0("G", 1:5))],
       xlab = "Orden cronologico de los archivos", ylab = "Media terminal (grados)",
       main = "Orden de adquisicion y grupos")
  boxplot(angulo_terminal_deg ~ grupo, data = ensayos, col = "grey90",
          xlab = "Grupo", ylab = "Media terminal (grados)", main = "Cinco ensayos por grupo")
  offsets <- seq(-0.15, 0.15, length.out = 5)
  for (i in 1:5) {
    points(i + offsets, ensayos$angulo_terminal_deg[ensayos$grupo == paste0("G", i)],
           pch = 19, col = colores[i])
  }
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  plot(resumen$densidad_g_cm3, resumen$media_deg, pch = 19,
       ylim = range(resumen[, c("ic95_inferior_deg", "ic95_superior_deg")]),
       xlab = "Densidad de referencia (g/cm3)", ylab = "Media terminal (grados)",
       main = "Media e IC95% entre ensayos")
  arrows(resumen$densidad_g_cm3, resumen$ic95_inferior_deg,
         resumen$densidad_g_cm3, resumen$ic95_superior_deg,
         angle = 90, code = 3, length = 0.05)
  plot(ensayos$angulo_terminal_deg, ensayos$densidad_g_cm3, pch = 19,
       col = colores[match(ensayos$grupo, paste0("G", 1:5))],
       xlab = "Media terminal (grados)", ylab = "Densidad de referencia (g/cm3)",
       main = "Calibracion lineal exploratoria")
  abline(modelo, lwd = 2)
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
  for (g in paste0("G", 1:5)) {
    x <- ensayos$angulo_terminal_deg[ensayos$grupo == g]
    qqnorm(x, main = paste("Q-Q", g, "(n=5)")); qqline(x)
  }
  plot.new()
  text(0.5, 0.5, "n=5 limita el diagnostico de normalidad.\nNo se demuestra independencia con un test.")
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  plot(fitted(modelo), residuals(modelo), xlab = "Densidad ajustada (g/cm3)",
       ylab = "Residuo (g/cm3)", main = "Residuos de calibracion", pch = 19)
  abline(h = 0, lty = 2)
  plot(seq_len(nrow(ensayos)), residuals(modelo), pch = 19,
       xlab = "Orden de ensayo", ylab = "Residuo (g/cm3)", main = "Residuos y orden")
  abline(h = 0, lty = 2)
  for (modo in c("ensayo", "densidad")) {
    z <- cv[cv$modo == modo, ]
    plot(z$densidad_referencia_g_cm3, z$densidad_predicha_g_cm3,
         xlab = "Referencia (g/cm3)", ylab = "Prediccion fuera de muestra (g/cm3)",
         main = paste("Excluir", modo), pch = 19,
         xlim = range(c(z$densidad_referencia_g_cm3, z$densidad_predicha_g_cm3)),
         ylim = range(c(z$densidad_referencia_g_cm3, z$densidad_predicha_g_cm3)))
    abline(0, 1, lty = 2)
  }
}
graficar()

fmt <- function(x) formatC(x, digits = 6, format = "g")
tabla_grupos <- c(
  "| Grupo | Densidad (g/cm3) | n | Media (grados) | SD entre ensayos | IC95% |",
  "|---|---:|---:|---:|---:|---|",
  vapply(seq_len(nrow(resumen)), function(i) {
    z <- resumen[i, ]
    sprintf("| %s | %.4f | %d | %.6f | %.6f | [%.6f; %.6f] |",
            z$grupo, z$densidad_g_cm3, z$n, z$media_deg, z$sd_entre_ensayos_deg,
            z$ic95_inferior_deg, z$ic95_superior_deg)
  }, character(1))
)
n_estables <- sum(ensayos$estabilidad_orientativa)
memoria <- c(
  "# Borrador de memoria - Modulo III",
  "",
  "## Titulo y estado",
  "Analisis estadistico de la repetibilidad y discriminacion de densidades mediante la orientacion de una capsula flotante.",
  "Borrador computacional, no entrega academica aprobada ni acreditacion del curso. No se han incorporado las autorizaciones del director ni de los coordinadores.",
  "",
  "## 1. Objetivo",
  "Cuantificar la variabilidad entre repeticiones y explorar si la orientacion terminal permite distinguir cinco condiciones de densidad. Evaluar una calibracion lineal sencilla sin asumir un modelo dinamico de primer orden.",
  "",
  "## 2. Datos y procedencia",
  paste("Fuente: BrewGold/Doctorado-densimetro, revision", revision),
  "Archivo: references/tfm/datos_udp.zip; 25 CSV. Asignacion explicita en datos/ensayos.csv, contrastada con Tabla C.1 del TFM; densidades de la Tabla 5.1.",
  "TFM: Estimacion de densidad de fluidos mediante orientacion angular, Fernando Gold Araoz, Rev04, 2026. Estos datos se reutilizan con atribucion; el presente trabajo es un reanalisis y no una nueva campana experimental.",
  "Cinco grupos de cinco ensayos; concentraciones de NaCl de 0, 2, 5, 12 y 26 g/L. Densidades de referencia tabuladas: 0.9973, 0.9988, 1.0011, 1.0064 y 1.0171 g/cm3; no se han verificado mediante una medicion independiente.",
  "Se analiza ang_deg en sus grados originales. No se reproduce la transformacion angular ni los parametros ajustados en radianes del TFM.",
  "",
  "## 3. Metodo",
  "Unidad experimental: el archivo/ensayo, segun el diseno documentado del TFM. Las miles de muestras temporales NO se tratan como repeticiones independientes.",
  "Se detiene el analisis ante columnas ausentes, valores no finitos o cobertura terminal insuficiente. Se separan segmentos ante diferencias temporales <= 0 y se utiliza solamente el ultimo segmento, sin ordenar artificialmente ni interpolar los datos.",
  "Medida principal: media aritmetica de ang_deg en (t_final - 10 s, t_final]. Sensibilidad: ventanas terminales de 5, 15 y 20 s. La media pondera muestras y supone que las perdidas no sesgan sistematicamente la orientacion.",
  "Cobertura minima para calcular: al menos 90% de muestras esperadas segun el intervalo mediano y ningun hueco mayor de 1 s en cada ventana. Huecos >100 ms activan una advertencia de cobertura; se mantienen los ensayos explicitamente senalados, sin relajar el indicador orientativo ni interpolar.",
  "Diagnostico orientativo de estabilidad: deriva lineal en 10 s <= 0.1 grados en valor absoluto, diferencia maxima entre medias de las cuatro ventanas <= 0.1 grados y cobertura sin huecos >100 ms. Son umbrales de trabajo, no especificaciones metrologicas validadas; ningun ensayo se excluye por este indicador.",
  "No se afirma que la media terminal sea el equilibrio fisico si los diagnosticos no lo apoyan. La SD temporal refleja fluctuacion intraensayo; no es la SD de repetibilidad ni se divide por la raiz del numero de muestras para inferencia.",
  "Por grupo: media, SD muestral entre cinco ensayos, SEM = SD/sqrt(5), IC95% = media +/- t(0.975,4)*SEM. Intervalos condicionados a independencia y distribucion razonablemente compatible con normalidad.",
  "Comparacion: ANOVA de Welch sin exigir igualdad de varianzas. Comparaciones t de Welch para los 10 pares, con p ajustadas por Holm. Los IC95% de las diferencias son nominales, no simultaneos. Nivel alfa 0.05.",
  "Diagnosticos: Q-Q y Shapiro-Wilk por grupo (baja potencia con n=5), Brown-Forsythe de varianzas y Kruskal-Wallis exploratorio. Este ultimo contrasta distribuciones/rangos; no debe interpretarse automaticamente como contraste de medianas cuando difieren las formas.",
  "Calibracion exploratoria: regresion OLS densidad ~ angulo terminal, con 25 valores. Validacion dejando fuera un ensayo y dejando fuera una densidad completa (cinco ensayos). La primera evalua nuevas repeticiones de condiciones conocidas; la segunda niveles no utilizados. Excluir G1/G5 implica extrapolacion en densidad de referencia.",
  "No se selecciona entre modelos mediante los mismos resultados de validacion ni se asume que R2 alto garantice capacidad predictiva. No se propaga aqui la incertidumbre de las densidades tabuladas.",
  "",
  "## 4. Resultados calculados",
  paste("Reinicios/duplicados temporales detectados:", sum(calidad$reinicios_o_duplicados),
        "; archivos afectados:", sum(calidad$reinicios_o_duplicados > 0), "."),
  paste("Muestras previas al segmento final no utilizadas:", sum(calidad$muestras_previas_no_usadas), "."),
  paste("Huecos >50 ms en registros originales:", sum(calidad$huecos_mayores_50ms), "."),
  paste("Ensayos que cumplen el diagnostico orientativo de estabilidad:", n_estables, "de 25."),
  paste("Hueco maximo entre las ventanas terminales evaluadas:",
        fmt(max(ensayos$max_hueco_ventanas_s)), "s."),
  paste("Ensayos con advertencia de cobertura:",
        paste(ensayos$archivo[!ensayos$cobertura_orientativa], collapse = ", ")),
  "Una advertencia por huecos no implica que la capsula sea fisicamente inestable.",
  if (n_estables < 25) "ADVERTENCIA: revisar los ensayos marcados FALSE; los contrastes describen medias terminales y no prueban diferencias entre equilibrios fisicos." else
    "Los 25 ensayos cumplen el criterio orientativo, lo que no demuestra por si solo equilibrio fisico.",
  "",
  tabla_grupos,
  "",
  paste("ANOVA de Welch: F =", fmt(welch$statistic), "; gl =", fmt(welch$parameter[1]),
        "y", fmt(welch$parameter[2]), "; p =", fmt(welch$p.value), "."),
  paste("Kruskal-Wallis exploratorio: H =", fmt(kw$statistic), "; p =", fmt(kw$p.value), "."),
  paste("Pares con p de Holm < 0.05:", sum(pares$p_holm < 0.05), "de 10.",
        "No confundir diferencias estadisticas con separacion fiable de medidas individuales."),
  paste("Pares sin evidencia suficiente tras Holm:",
        if (all(pares$p_holm < 0.05)) "ninguno." else
          paste(paste(pares$grupo_a[pares$p_holm >= 0.05],
                      pares$grupo_b[pares$p_holm >= 0.05], sep = "-"), collapse = ", ")),
  paste("Diagnostico de normalidad, p de Shapiro-Wilk por grupo (sin ajuste, exploratorias):",
        paste(paste(normalidad$grupo, fmt(normalidad$p), sep = "="), collapse = "; "), "."),
  paste("Diagnostico de varianzas Brown-Forsythe: p =", fmt(bf$`Pr(>F)`[1]), "."),
  "Los diagnosticos deben leerse junto con los Q-Q, la elevada dispersion de G1 y el tamano muestral; no validan por si solos la inferencia parametrica. No se eliminan valores por parecer atipicos.",
  paste("Calibracion: densidad [g/cm3] =", fmt(coef(modelo)[1]), "+",
        fmt(coef(modelo)[2]), "* angulo_terminal [grados]. R2 de ajuste =",
        fmt(summary(modelo)$r.squared), "."),
  vapply(seq_len(nrow(metricas)), function(i) {
    z <- metricas[i, ]
    paste("Validacion excluyendo", z$modo, ": RMSE =", fmt(z$rmse_g_cm3),
          "; MAE =", fmt(z$mae_g_cm3), "; sesgo =", fmt(z$sesgo_g_cm3), "g/cm3.")
  }, character(1)),
  "Las p de los coeficientes OLS se exportan como salida convencional; no son evidencia confirmatoria de una ley fisica ni sustituyen la validacion.",
  "",
  "## 5. Discusion y limitaciones",
  "La campana agrupa densidades en orden cronologico. Densidad y orden de adquisicion quedan confundidos: no puede atribuirse causalmente toda diferencia al fluido ni descartarse deriva instrumental con estos datos.",
  "La variable temp corresponde al sensor; no debe usarse como temperatura medida del fluido. El TFM indica temperatura ambiente aproximada de 24 grados C, pero los valores del sensor no verifican esa condicion.",
  "La independencia entre ensayos procede del protocolo declarado y requiere revisar cuaderno de laboratorio, preparacion de disoluciones y reinicios; no la demuestra una prueba estadistica. Si se reutilizo una misma disolucion por grupo, las cinco repeticiones no son cinco preparaciones independientes.",
  "Con n=5 por grupo, la inferencia es limitada. Un p grande no demuestra igualdad ni ausencia de relevancia practica. Los IC no incluyen incertidumbre de referencia, deriva ni sesgos sistematicos.",
  "El criterio terminal evita imponer una relajacion monotona de primer orden a transitorios oscilatorios. No identifica frecuencias, amortiguamiento ni viscosidad; estos objetivos quedan fuera del alcance.",
  "La validacion interna no sustituye una campana independiente con orden aleatorizado, temperatura de fluido registrada y densidad de referencia medida. No se convierte el RMSE en resolucion o incertidumbre metrologica certificada.",
  "",
  "## 6. Conclusiones provisionales",
  if (welch$p.value < 0.05)
    "Los datos muestran evidencia de diferencias entre las medias terminales de al menos dos condiciones bajo los supuestos indicados; no prueban que todas las densidades sean distinguibles ni una relacion causal exclusiva con densidad." else
    "El contraste global no aporta evidencia suficiente de diferencias de medias terminales al nivel 0.05; esto no demuestra equivalencia entre condiciones.",
  paste("La estabilidad terminal debe interpretarse a partir de los", n_estables,
        "ensayos que cumplen el criterio orientativo y del analisis de sensibilidad, no asumirse."),
  "La calibracion queda como exploratoria: su utilidad depende del error fuera de muestra y de la tolerancia requerida, que aun no se ha definido. Debe revisarse con el director antes de formular conclusiones de prestaciones.",
  "",
  "## 7. Reproducibilidad y anexos",
  "Ejecutar Rscript tareas/modulo_III/scripts/analizar.R desde la raiz del repositorio. Solo requiere R base; no instala paquetes. Ver README.md para trabajo sin conexion y pruebas.",
  paste("MD5 del ZIP para integridad reproducible (no firma de seguridad):", md5_esperado),
  "Anexos: scripts/analizar.R y scripts/funciones.R contienen todo el codigo; resultados/*.csv contiene trazabilidad, estadisticos, diagnosticos y predicciones; resultados/figuras.pdf contiene series y graficos estadisticos; resultados/sessionInfo.txt registra el entorno.",
  "",
  "## 8. Bibliografia y fuentes",
  "UNED (2024-2025). Guia de estudio: Gestion y Analisis de Datos Cientificos. Modulo III de las actividades formativas transversales. Archivo ../Guia_Modulo_III.pdf.",
  "Gold Araoz, F. (2026). Estimacion de densidad de fluidos mediante orientacion angular. TFM, Rev04. Tablas 5.1 y C.1. Archivo references/tfm/Densidad_orientacion_angular-F_Gold_Rev04.pdf del repositorio fuente.",
  "BrewGold. Doctorado-densimetro. Datos experimentales y docs/PRELIMINARES_INVESTIGACION.md. Consultar revision fijada para trazabilidad, no confundir datos historicos con simulaciones.",
  "R Core Team. R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing. Version exacta en sessionInfo.txt. https://www.R-project.org/",
  "Documentacion de R, paquete stats: oneway.test, t.test, p.adjust, shapiro.test, kruskal.test y lm.",
  "",
  "## 9. Pendientes academicos",
  "Completar fechas de realizacion y revisar interpretacion, figuras y atribucion con el director. Confirmar con los coordinadores si procede una entrega o solo documentar el curso ya realizado; no presentar este reanalisis como trabajo previamente aprobado.",
  "La guia disponible corresponde a 2024-2025 e indica 15 de mayo: confirmar convocatoria y plazo aplicables, sin extrapolar esa fecha automaticamente al curso actual."
)
memoria <- c(memoria, "", "## Anexo A. Codigo R completo", "",
             "Los archivos originales adjuntos son la version ejecutable; el PDF presenta el codigo como texto de consulta.")
for (nombre in c("funciones.R", "analizar.R")) {
  memoria <- c(memoria, "", paste("### scripts/", nombre, sep = ""), "", "```r",
               readLines(file.path(raiz, "scripts", nombre), encoding = "UTF-8"), "```")
}
writeLines(memoria, file.path(salidas, "memoria_borrador.md"), useBytes = TRUE)
pdf_texto(memoria, file.path(salidas, "memoria_borrador.pdf"))
propuesta <- readLines(file.path(raiz, "propuesta.md"), encoding = "UTF-8")
pdf_texto(propuesta, file.path(salidas, "propuesta_borrador.pdf"))
writeLines(trimws(capture.output(sessionInfo()), which = "right"),
           file.path(salidas, "sessionInfo.txt"))
writeLines(c(paste("revision_fuente:", revision), paste("url_zip:", url),
             paste("md5_zip:", md5_esperado), "ventana_principal_s: 10",
             "ventanas_sensibilidad_s: 5, 15, 20",
             "tolerancia_deriva_10s_deg: 0.1", "tolerancia_ventanas_deg: 0.1",
             "max_hueco_advertencia_s: 0.1", "max_hueco_calculo_s: 1",
             "fraccion_minima_muestras_esperadas: 0.9", "alfa: 0.05"),
           file.path(salidas, "procedencia.txt"))
message("Resultados guardados en ", salidas)
if (n_estables < 25) message("ATENCION: ", 25 - n_estables,
                            " ensayos con advertencias de cobertura/estabilidad.")
