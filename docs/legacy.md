# Interfaz anterior con VCF separados

Imputa variantes llamadas en un panel de referencia usando genotipos SNP y
evalúa su precisión por alelo mediante validación cruzada anidada.
Acepta nombres de cromosomas de cualquier especie **diploide**. Versión inicial
0.1.0: ridge con dosages compatibles por locus; no incluye el descubrimiento
de variantes, genotipado desde lecturas ni los modelos de predicción BLUPF90.

## Entradas

- `--snps-vcf`: SNPs bialélicos con `GT`, en referencias y animales a imputar.
- `--targets-vcf`: SNPs o SVs con `GT` en los animales de referencia.
- `--predictor-pool`: lista opcional de IDs o claves `CHROM:POS:REF:ALT`.
  Define el panel disponible; en el artículo corresponde a los 18.007 SNPs.
  Sin lista, se usa el panel completo del VCF de SNPs.

Los IDs de muestras deben coincidir exactamente. Todos los animales del VCF
de variantes deben estar en el de SNPs. Los animales adicionales del VCF SNP
son los destinatarios de la imputación. Ambos VCF deben usar la misma referencia
y los mismos nombres de cromosoma. Se admiten VCF de texto o comprimidos que
pueda leer pysam; no hace falta un índice para la lectura secuencial.

## Instalación y uso

```bash
python -m pip install .

# Comprobar entradas y representación sin ajustar modelos
sv-impute --snps-vcf SNPs.vcf.gz --targets-vcf SVs.vcf.gz \
  --predictor-pool panel18K.txt --check-only --out runs/check

# Evaluar cada variante del cromosoma 23; diez particiones independientes
sv-impute --snps-vcf SNPs.vcf.gz --targets-vcf SVs.vcf.gz \
  --predictor-pool panel18K.txt --chromosome 23 \
  --predictor-mode nearest --neighbors 64 \
  --outer-folds 5 --inner-folds 3 --repeats 10 --out runs/chr23
```

`--chromosome` puede repetirse; si se omite se evalúan todos los cromosomas
del VCF objetivo. `--target-ids` permite elegir variantes concretas. Cada
directorio de salida debe ser nuevo o estar vacío.

| Modo | Predictores por variante |
| --- | --- |
| `nearest` | Hasta `--neighbors` SNPs cercanos del mismo cromosoma, dentro del panel |
| `chromosome` | Todos los SNPs del panel en ese cromosoma |
| `pool` | Todos los SNPs del panel definido por el usuario |

El tamaño del vecindario es una decisión del modelo, no un filtro de calidad.
Si el panel incluye la posición objetivo entre los predictores elegidos, el
programa registra el conflicto y no ajusta ese locus: debe resolverse el panel
para evitar predecir una variante a partir de sí misma.

## Método y resultados

Ridge selecciona una penalización por alelo en los folds internos. La media
para SNPs ausentes y la estandarización se estiman exclusivamente en el
entrenamiento de cada fold. Los folds externos proporcionan predicciones
fuera de muestra. Se guardan los splits, parámetros, coeficientes y hashes.
Las predicciones se acotan a [0,2] y se proyectan conjuntamente para que la suma
de copias ALT no supere dos. Se conservan resultados antes y después de la
proyección. La proyección no resuelve errores del genotipo de referencia.

- `accuracy.tsv`: Pearson **con signo** entre GT 0/1/2 y dosage predicho,
  MSE y número de pares por alelo, repetición y tipo de predicción.
- `accuracy_summary.tsv`: media y desviación de r entre repeticiones válidas.
  Las repeticiones comparten animales; la desviación no es un intervalo de confianza.
- `status.tsv`: cada locus solicitado, incluidos los que no se pudieron calcular.
- `predictors.tsv`: SNPs usados para cada locus.
- `imputed_repeatNN.vcf`: DS por ALT en animales nuevos; GT permanece ausente.
- `loci/`: llamadas originales, estados por muestra y resultados completos.
- `splits.json`, `run.json`: particiones, opciones y procedencia.

No se aplican cortes de MAF, calidad, LD, precisión ni significancia. Las
correlaciones de genotipos constantes quedan indefinidas, con su motivo.
Esto es precisión de **dosage**, no precisión de genotipos discretizados.

Un VCF multialélico conserva un locus y todos sus ALT. Los registros separados
con la misma posición se revisan conjuntamente: si sus copias ALT exceden dos,
no se ajusta el locus ni se elige un ALT arbitrario. Si REF/END difieren o un ALT
está duplicado, se conserva la evidencia y se informa que falta resolver su
representación; ese locus no se reconstruye en el VCF de salida. La compatibilidad
aritmética por sí sola no demuestra que las llamadas sean biológicamente correctas.

La versión inicial procesa cada locus secuencialmente y mantiene el panel SNP
en memoria. Es adecuada para paneles definidos por el usuario; no está optimizada
para cargar millones de SNPs simultáneamente. No admite genotipos poliploides.

## Particiones explícitas

`--splits-json` acepta el formato de `splits.json`: `sample_ids` en el orden
exacto del VCF objetivo y `repeats`, una lista de particiones. Cada repetición
contiene, por fold externo, los folds internos seguidos del externo. Cada
partición tiene índices `fit` y `evaluate`, empezando en cero. Se comprueba que
no haya solapamiento y que cada muestra se evalúe una vez por repetición.

Los folds aleatorios por individuo son el valor inicial. Para medir transferencia
entre familias o poblaciones se deben proporcionar particiones apropiadas.
La división automática no identifica parentescos.

## Pruebas

```bash
python -m unittest discover -s tests -v
```

Las pruebas usan datos sintéticos. No se incluyen genotipos, fenotipos,
credenciales ni rutas privadas del artículo en este repositorio.
