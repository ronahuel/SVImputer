# SVImputer

SVImputer predicts structural variant dosages from a mixed SNP/SV reference VCF and a low-density SNP panel in VCF or PLINK BIM format. It compares nearby SNP counts by cross-validation and plots mean imputation accuracy. Each run also describes the reference VCF. To generate only the descriptive figures, run:

```bash
python -m pip install .
sv-impute describe --reference examples/allelic_series.vcf --out example_profile
```

The example includes biallelic, triallelic and four-allele loci. One triallelic site is split over two VCF lines. The profile reports the alleles listed in the VCF and accepts sites-only VCFs.

See the example plots for [SV types](examples/figures/sv_types.png), [SV lengths](examples/figures/sv_lengths.png) and [alleles per locus](examples/figures/alleles_per_locus.png).

## Guía de uso

Predice el dosage de cada alelo ALT a partir de los SNPs más cercanos del panel LD.
Con un VCF de referencia y un panel LD en VCF o BIM obtienes el recuento de marcadores por cromosoma, las predicciones y
una figura que muestra con cuántos SNPs vecinos se alcanza la mayor precisión media.

## Probar con tus datos

En Linux o WSL, desde la carpeta del repositorio, con Python 3.9 o superior:

```bash
python -m pip install .
sv-impute --reference referencia.vcf.gz --ld-panel-vcf panel_LD.vcf.gz --out resultados_vcf
sv-impute --reference referencia.vcf.gz --ld-panel-bim panel_LD.bim --out resultados_bim
```

- **Referencia.** VCF ordenado por cromosoma y posición, con SNPs y SVs juntos y genotipos diploides en `FORMAT/GT`. Se aceptan llamadas faseadas o sin fase.
- **`--ld-panel-vcf`** acepta un panel de SNPs en VCF con solo sitios o con genotipos. Si contiene animales adicionales con GT, también se predicen sus dosages. `--ld` sigue disponible como alias.
- **`--ld-panel-bim`** acepta un mapa PLINK de texto (`.bim` o `.bim.gz`). Para validar con los animales del VCF de referencia basta el BIM. El mapa identifica los SNPs y el programa toma sus genotipos del VCF de referencia. Elige una sola opción de panel por corrida.
- Las columnas del mapa son `CHR ID cM BP A1 A2` (también se admite omitir cM). El emparejamiento usa el cromosoma, la posición física BP y ambos alelos en cualquier orden. Un alelo ausente (`0` o `.`) permite usar la posición y el alelo conocido si identifican un único SNP. El informe `ld_matches.tsv` conserva la fila y el motivo si hay discrepancia o ambigüedad. [Formato BIM de PLINK](https://www.cog-genomics.org/plink/1.9/formats#bim)
- Usa el mismo ensamblaje y los mismos nombres de cromosoma en ambos archivos. Se admite VCF de texto o comprimido con **bgzip**. La lectura es secuencial, por lo que el índice es opcional.

Para comparar ambos formatos con el panel sintético incluido, ejecuta:

```bash
sv-impute --reference examples/reference.vcf --ld-panel-vcf examples/ld.vcf --neighbors 1 2 4 --out demo_vcf
sv-impute --reference examples/reference.vcf --ld-panel-bim examples/ld.bim --neighbors 1 2 4 --out demo_bim
```

Para cambiar los tamaños de vecindario, limitar el análisis a un cromosoma o repetir la validación:

```bash
sv-impute --reference referencia.vcf.gz --ld-panel-vcf panel_LD.vcf.gz \
  --chromosome chr1 --neighbors 4 8 16 32 64 --repeats 3 --out resultados_chr1
```

## Cómo leer los resultados

Abre **`accuracy_by_neighbors.png`**. El eje horizontal muestra el máximo de SNPs vecinos y el vertical la correlación de Pearson **con signo** entre dosage predicho fuera de muestra y copias ALT observadas (0, 1 o 2). Los puntos destacados son los máximos entre los tamaños probados. Si dos tamaños empatan, aparecen ambos. El máximo resume esta comparación exploratoria en la referencia. Una prueba en animales independientes permite medir su rendimiento en otra población.

| Archivo | Contenido |
| --- | --- |
| `reference_profile/summary.json` y `reference_profile/chromosomes.tsv` | Versión del VCF, muestras, cromosomas con registros, filas por categoría y loci por posición |
| `reference_profile/sv_types.png` y `sv_types.tsv` | Número de entradas ALT de SV por tipo y cuántas tienen longitud conocida |
| `reference_profile/sv_lengths.png` y `sv_alleles.tsv.gz` | Distribución acumulada de la longitud en pb por tipo y detalle por entrada, con el método o el motivo de ausencia |
| `reference_profile/alleles_per_locus.png` y `alleles_per_locus.tsv` | Loci bialélicos, trialélicos, con cuatro o más alelos y sin resolver. Incluye el conjunto completo y el subconjunto de loci con SVs |
| `chromosomes.tsv` | SNPs, SVs, indels y otros registros de referencia, SNPs LD recibidos y emparejados, loci evaluados por cromosoma |
| `window_accuracy.tsv` | Precisión media global y por cromosoma, número de loci comparables, rango real de SNPs utilizados y mejores tamaños |
| `accuracy.tsv.gz` | Pearson, número de pares o motivo de cálculo indefinido por locus, ALT, tamaño y repetición |
| `dosage_cv.tsv.gz` | Dosages de validación para los individuos de referencia. Cada columna es un individuo |
| `dosage_ld.tsv.gz` | Dosages de individuos nuevos, si están en el VCF LD. Para el tamaño preferido se selecciona su valor en la columna `k` |
| `ld_matches.tsv`, `neighbors.tsv.gz` | Estado de cada SNP suministrado y posiciones vecinas por locus en orden de distancia |
| `splits.json`, `run.json` | Particiones, parámetros, número de muestras y tamaños con mejor media |

**Selección de vecinos.** Se prueban por defecto 1, 2, 4, 8, 16, 32, 64 y 128 SNPs del mismo cromosoma. Si un cromosoma tiene menos SNPs se usan los disponibles. Los tamaños por encima del máximo del panel quedan representados por ese máximo. En `neighbors.tsv.gz`, los primeros `k` SNPs son los de esa ventana. A igual distancia se toma primero el de menor posición. Un SNP en la posición objetivo nunca sirve de predictor para ese locus.

**Variantes.** Por defecto se evalúan las SVs (anotadas con `SVTYPE`, simbólicas o breakends) y los indels de secuencia de cualquier longitud. Usa `--targets sv` para las SVs identificables por anotación/representación, `--targets snp` para los SNPs ajenos al panel o `--targets all` para todas las variantes objetivo. Los SNPs LD se emparejan por cromosoma, posición y alelos. Si REF y ALT están invertidos, se reorientan los dosages del panel. Las ausencias, discrepancias y ambigüedades quedan en `ld_matches.tsv`.

**Validación.** Todos los tamaños comparten las mismas particiones. Ridge selecciona la penalización por alelo en los folds internos (3 por defecto). Los externos (5 por defecto) generan una predicción fuera de muestra por animal en cada repetición. La media y la escala de los SNPs se aprenden en el entrenamiento. En los loci multialélicos, la suma de dosages ALT se limita a dos copias. Las llamadas incompatibles quedan identificadas. Si un animal carece de todos sus SNPs locales, su dosage se registra como ausente. No se aplican filtros de MAF, calidad, LD ni precisión.

**Promedio.** Para cada tamaño se promedia primero entre alelos y repeticiones dentro de cada locus, y luego entre loci con igual peso. La figura usa los mismos pares alelo–repetición calculables en todos los tamaños de cada locus (`mean_r_paired`). La tabla también entrega la media de todos los valores disponibles para cada tamaño (`mean_r_available`). Un Pearson indefinido queda en blanco con su motivo. Si ningún locus permite la comparación, el programa lo indica en la figura.

`--chromosome` puede repetirse. Usa `--splits-json` para particiones propias (p. ej., por familia), en el formato de `splits.json`. Las divisiones automáticas son aleatorias por individuo. El programa mantiene el panel LD en memoria y lee las variantes objetivo por locus. El tiempo y el tamaño de la salida crecen con los loci, vecindarios y repeticiones.

### Cómo se describe el VCF

El informe abarca el VCF de referencia dentro de los cromosomas solicitados, incluidos los marcadores ajenos al panel LD. En el comando independiente `describe`, los archivos se guardan directamente en `--out`.

Para contar alelos, cada combinación de cromosoma y posición representa un locus. Con una REF común se suma esa REF y cada ALT distinta, aunque aparezca en varias líneas. Los marcadores de sitio ausente, alelo solapante o alelo no especificado (`.`, `*`, `<*>`, `<NON_REF>`) no añaden alelos. Las REF discordantes, un ALT simbólico repetido con extensiones o longitudes distintas y los sitios sin ALT explícita quedan en `unresolved_loci.tsv.gz`. El conteo describe el catálogo de alelos representados, independientemente de su frecuencia en los genotipos.

Los gráficos de tipo y longitud cuentan cada entrada ALT estructural del VCF. Se usa el tipo simbólico, la notación de breakend o la etiqueta SVTYPE. Los indels literales sin anotación se describen como indels en el resumen. Para la longitud se usa el valor absoluto de SVLEN por ALT, la diferencia entre las longitudes de secuencia en INS/DEL o, en DEL/DUP/INV/CNV simbólicas de un solo ALT, el intervalo de referencia cuando está definido. Las longitudes desconocidas y los breakends siguen presentes en los conteos y en la tabla. [Especificación VCF](https://samtools.github.io/hts-specs/VCFv4.5.pdf)

La [interfaz anterior](docs/legacy.md) sigue disponible. Para ejecutar las pruebas, usa `python -m unittest discover -s tests -v`.
