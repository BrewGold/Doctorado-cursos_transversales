# Propuesta de trabajo - Modulo III

## Titulo

Analisis estadistico de la repetibilidad y discriminacion de densidades
mediante la orientacion de una capsula flotante.

## Estado y contexto

Borrador para revisar con el director de tesis. El curso ya se realizo
anteriormente; este documento organiza un posible reanalisis de datos existentes
y no afirma que sea la propuesta aprobada en aquella ocasion. Cualquier uso como
nueva entrega requiere confirmar su pertinencia y obtener las autorizaciones
academicas correspondientes.

El trabajo se vincula a la investigacion sobre un densimetro basado en la
orientacion de una capsula flotante con una unidad de medida inercial. Se dispone
de una campana historica del TFM con 25 registros temporales: cinco ensayos en cada
una de cinco condiciones. Los datos se conservan en Doctorado-densimetro y se
utilizaran sin alterar los originales.

## Objetivo del estudio

El objetivo principal es cuantificar la variabilidad entre repeticiones y
explorar si las diferencias de orientacion terminal entre condiciones superan
esa variabilidad. Como objetivo secundario, se evaluara una calibracion lineal
de densidad a partir de la orientacion, distinguiendo bondad de ajuste y
capacidad predictiva fuera de muestra.

Se estudiaran las referencias 0.9973, 0.9988, 1.0011, 1.0064 y 1.0171 g/cm3,
correspondientes a concentraciones de NaCl de 0, 2, 5, 12 y 26 g/L segun el
TFM. Son valores tabulados, no densidades verificadas mediante una nueva
medicion. No se pretende identificar viscosidad ni resolver la dinamica
completa de la capsula.

## Metodo de desarrollo

1. Trazabilidad y control de calidad. Fijar una revision del repositorio fuente,
verificar el archivo ZIP y registrar la correspondencia explicita de cada CSV
con su grupo y repeticion mediante las tablas del TFM. Comprobar columnas,
valores finitos, continuidad temporal y cobertura de adquisicion. Ante un
reinicio temporal, separar segmentos y analizar solo el ultimo; registrar las
muestras anteriores no utilizadas sin reordenar ni interpolar.

2. Extraccion por ensayo. Calcular la media de ang_deg en los ultimos 10 s
como variable terminal, conservando su convencion original en grados.
Comprobar deriva lineal y sensibilidad a ventanas de 5, 15 y 20 s. Aplicar
indicadores orientativos de cobertura y estabilidad, declarando sus umbrales,
sin excluir automaticamente valores que no los cumplan. No denominar equilibrio
fisico a la media terminal sin respaldo de los diagnosticos.

3. Analisis descriptivo. Representar las series, los cinco valores individuales
por grupo, medias, desviaciones estandar e intervalos de confianza del 95%.
Calcular la incertidumbre de la media usando cinco ensayos, no el numero de
muestras temporales. Distinguir ruido dentro de un registro y repetibilidad
entre registros.

4. Comparacion de grupos. Emplear ANOVA de Welch y comparaciones t de Welch
con ajuste de Holm para los diez pares. Inspeccionar Q-Q, normalidad y varianzas,
reconociendo la limitada informacion de cinco observaciones por grupo. Incluir
Kruskal-Wallis como sensibilidad sobre distribuciones, sin interpretarlo
automaticamente como una comparacion de medianas.

5. Calibracion y validacion. Ajustar un modelo lineal exploratorio densidad
en funcion del angulo terminal, revisar residuos y estimar RMSE, MAE y sesgo.
Validar dejando fuera ensayos completos y, por separado, grupos de densidad
completos. Explicar que excluir una densidad extrema introduce extrapolacion.
No convertir errores de prediccion en prestaciones metrologicas certificadas.

## Fases y fechas de realizacion

El calendario siguiente es relativo y propuesto; las fechas naturales y el
inicio se acordaran con el director. No presenta actividades como aprobadas
ni retrospectivamente realizadas en fechas inventadas.

Fase 1, dia 1, 1 hora: revisar procedencia, manifiesto y calidad de datos.
Fase 2, dia 2, 2 horas: programar extraccion terminal y analisis descriptivo.
Fase 3, dia 3, 2 horas: comparaciones, diagnosticos y calibracion.
Fase 4, dia 4, 2 horas: interpretar resultados y preparar memoria y anexos.

Son siete horas orientativas de trabajo practico, en consonancia con la
planificacion de la guia; revisar la estimacion si la depuracion o la revision
requieren mas tiempo. El estudio del texto base y las revisiones academicas
son adicionales. La guia disponible, de 2024-2025, indica 15 de mayo como plazo;
no se asume que esa fecha sea aplicable a una convocatoria actual.

## Medios y justificacion

Se empleara R base por su soporte de analisis estadistico, modelos lineales,
graficos y automatizacion reproducible. No se requieren paquetes CRAN
adicionales ni adquisicion de nuevas medidas para este alcance. El codigo,
el manifiesto, los parametros, las tablas y los documentos se conservaran en
el repositorio de cursos transversales; los datos originales seguiran en su
repositorio de origen con referencia a una revision fija.

Los medios son pertinentes porque permiten reconstruir el flujo completo
desde registros experimentales hasta conclusiones, estudiar distribuciones,
comparar grupos y evaluar un modelo sin depender de resultados copiados
del TFM.

## Resultados previstos y limitaciones

Se entregaran codigo R, tablas CSV, un anexo de figuras PDF y una memoria PDF
con objetivos, analisis, resultados, conclusiones y bibliografia. Se informaran
tanto diferencias detectadas como comparaciones inconcluyentes y advertencias
de calidad, sin seleccionar solo resultados favorables.

La inferencia depende de la independencia de los ensayos y de sus supuestos
distribucionales. El orden de adquisicion esta confundido con la densidad.
La temperatura del sensor no equivale a temperatura del fluido. La reutilizacion
de disoluciones y la incertidumbre de los valores tabulados deben revisarse;
los registros por si solos no resuelven esos limites.

La aportacion frente al TFM es un reanalisis estadistico explicito de medias
terminales, calidad temporal y validacion predictiva, sin imponer el modelo
dinamico de primer orden ni atribuir de antemano todos los cambios a densidad.

## Referencias y revision

UNED (2024-2025). Guia de estudio: Gestion y Analisis de Datos Cientificos,
Modulo III de las actividades formativas transversales.

Gold Araoz, F. (2026). Estimacion de densidad de fluidos mediante orientacion
angular. TFM, Rev04. Tablas 5.1 y C.1.

R Core Team. R: A Language and Environment for Statistical Computing.
Documentacion de stats y version del entorno en resultados/sessionInfo.txt.

Antes de una eventual entrega: revision y autorizacion del director, aprobacion
de la propuesta por los coordinadores y revision final de la memoria por el
director, conforme al procedimiento de la guia.
