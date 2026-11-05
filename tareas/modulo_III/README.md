# Módulo III: Gestión y Análisis de Datos Científicos

Reanálisis reproducible en **R base** de 25 ensayos históricos del densímetro:
cinco repeticiones en cada una de cinco densidades de referencia. Se sigue el
alcance de la [guía del módulo](../Guia_Modulo_III.pdf): análisis de datos de
investigación, comparación de grupos, modelado y memoria con código y bibliografía.

**Estado:** borrador para revisión. No acredita aprobación académica ni sustituye
el trabajo realizado anteriormente. La reutilización de datos y la diferencia
respecto del TFM deben declararse. La guía es de 2024–2025; confirmar plazos y
requisitos vigentes con el director y los coordinadores si se pretende entregar.

## Ejecutar

Desde la raíz del repositorio:

```sh
Rscript tareas/modulo_III/scripts/analizar.R
Rscript tareas/modulo_III/tests/pruebas.R
```

Requiere R >= 4.1 con los paquetes estándar `stats`, `utils`, `tools`, `graphics`
y `grDevices`; no utiliza paquetes CRAN adicionales. El entorno exacto de la
ejecución publicada está en `resultados/sessionInfo.txt`.

La primera ejecución descarga el ZIP experimental de una revisión fija de
`BrewGold/Doctorado-densimetro` y verifica su MD5. Las siguientes utilizan esa
caché; no consultan `main`. El MD5 comprueba integridad reproducible, no autenticidad
criptográfica. La descarga puede fallar explícitamente si no hay acceso a GitHub.
Para trabajar sin conexión, guardar el ZIP de la URL indicada en
`resultados/procedencia.txt` como `datos/cache/datos_udp.zip`. No modificarlo.
Los CSV se extraen de nuevo del ZIP verificado en cada ejecución.

La caché está excluida de Git para no duplicar los datos originales. Las tablas,
los PDF y la memoria generada sí se versionan. El script sobrescribe únicamente
sus salidas conocidas; no modifica el repositorio doctoral ni sus datos.
`memoria_borrador.md` y los PDF son generados: cambiar el generador para conservar
ediciones de la memoria; editar `propuesta.md` para cambiar la propuesta.

## Archivos y resultados

| Ruta | Contenido |
|---|---|
| `datos/ensayos.csv` | Correspondencia explícita archivo–grupo–densidad–repetición |
| `propuesta.md` | Objetivos, método, fases, calendario relativo y medios |
| `scripts/analizar.R` | Descarga, análisis, figuras y generación de memoria |
| `scripts/funciones.R` | Preparación temporal, resúmenes, contrastes y validación |
| `tests/pruebas.R` | Casos sintéticos y comprobación de resultados publicados |
| `resultados/calidad.csv` | Reinicios, segmentos, muestras no utilizadas y huecos |
| `resultados/ensayos.csv` | Una fila por ensayo, ventanas y diagnósticos terminales |
| `resultados/resumen_grupos.csv` | Media, SD, SEM e IC95% entre cinco ensayos |
| `resultados/contrastes_globales.csv` | Welch y Kruskal–Wallis exploratorio |
| `resultados/comparaciones_pares.csv` | Diez comparaciones Welch y ajuste de Holm |
| `resultados/normalidad_por_grupo.csv` | Shapiro–Wilk exploratorio |
| `resultados/homogeneidad_varianzas.csv` | Brown–Forsythe |
| `resultados/sensibilidad_ventanas.csv` | Comparación de ventanas de 5, 10, 15 y 20 s |
| `resultados/calibracion_*.csv` | Coeficientes, predicciones excluidas y errores |
| `resultados/figuras.pdf` | Series, dispersión, Q-Q, residuos y calibración |
| `resultados/memoria_borrador.md` | Memoria con resultados, bibliografía y código completo |
| `resultados/memoria_borrador.pdf` | Versión PDF de lectura, con las tablas como texto |
| `resultados/propuesta_borrador.pdf` | Propuesta en PDF, pendiente de revisión |
| `resultados/procedencia.txt` | Revisión, URL, checksum y parámetros del análisis |

Las figuras son un anexo separado de la memoria PDF. Para una entrega, adjuntar
ambos documentos y el código, o maquetarlos conjuntamente tras la revisión.

## Decisiones metodológicas

- **Unidad experimental:** ensayo, no muestra temporal. La independencia se asume
  según el protocolo del TFM; no se demuestra con los CSV.
- **Tiempos:** separar ante `diff(t_ms) <= 0` y usar el último segmento. No ordenar
  muestras para ocultar reinicios, unir segmentos ni interpolar huecos.
- **Variable principal:** media de `ang_deg` en los últimos 10 s del segmento
  final. Se conserva la convención angular original en grados, sin equipararla
  a los parámetros en radianes del TFM.
- **Cobertura:** se detiene el cálculo si una ventana tiene menos del 90 % de las
  muestras esperadas por su intervalo mediano, menos de 10 muestras o un hueco
  mayor de 1 s. Se advierte si un hueco supera 100 ms. También se exige un segmento
  final de al menos 30 s. Se rechazan valores ausentes/no finitos en las columnas
  requeridas, sin imputación ni eliminación silenciosa.
- **Estabilidad orientativa:** deriva lineal en 10 s y variación de medias entre
  ventanas ≤ 0,1° cada una, y cobertura sin huecos > 100 ms. Estos umbrales son
  decisiones exploratorias, no tolerancias certificadas del instrumento.
  Los indicadores no excluyen ensayos automáticamente ni prueban equilibrio físico.
- **Incertidumbre:** SD y SEM entre ensayos; IC95% con t de Student y cuatro
  grados de libertad. No utilizar el tamaño de la serie temporal como `n`.
- **Comparación:** Welch global y diez pruebas t de Welch con p ajustadas por
  Holm. Los IC95% de diferencias son nominales, no simultáneos. Los diagnósticos
  de normalidad son limitados con cinco observaciones; Kruskal–Wallis es una
  sensibilidad sobre distribuciones, no una prueba automática de medianas.
- **Calibración:** modelo exploratorio `densidad ~ angulo_terminal`, sin selección
  oportunista de modelos. Excluir un ensayo evalúa condiciones conocidas; excluir
  una densidad evalúa niveles no utilizados y puede implicar extrapolación.

### Limitaciones que no deben omitirse

Las densidades son referencias tabuladas del TFM, no medidas independientes
verificadas aquí. El orden de adquisición está confundido con la densidad y la
variable `temp` es del sensor, no una medición del fluido. Deben revisarse
preparación de disoluciones, independencia y posibles derivas antes de interpretar
causalmente las diferencias. Los errores predictivos no son resolución ni
incertidumbre metrológica certificada.

## Procedencia

Fuente fijada: `BrewGold/Doctorado-densimetro`,
revisión `bcc97816212b48ea0aed2849c85d9b0065c223ed`.
Datos: `references/tfm/datos_udp.zip`.

La asignación de los 25 archivos procede de la Tabla C.1 y las densidades de la
Tabla 5.1 de *Estimación de densidad de fluidos mediante orientación angular*,
Fernando Gold Araoz, Rev04, 2026:
`references/tfm/Densidad_orientacion_angular-F_Gold_Rev04.pdf`.

El TFM se conserva como referencia histórica. Este trabajo no reutiliza sin
validación sus scripts ni asume que el modelo de primer orden sea válido para
los transitorios; consultar también `docs/PRELIMINARES_INVESTIGACION.md` en la
misma revisión del repositorio fuente.
