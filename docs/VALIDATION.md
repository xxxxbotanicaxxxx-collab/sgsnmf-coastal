# Validación de SGSNMF Coastal 0.1.0

Ejecutada el 4 de septiembre de 2026 en Windows, Python 3.12 y GNU Octave
10.3.0. Dependencias exactas en `requirements.lock.txt`. Los resultados
describen esta ejecución; no implican igualdad entre versiones de BLAS,
SciPy u Octave.

## Comparación con el original

El demo original completo (221 bandas, 100×100 píxeles, 9 endmembers,
20 dB) se procesó durante 100 iteraciones en ambos lenguajes, desactivando
el límite CPU para esta prueba. Python recibió exactamente X, W0, H0 y
la segmentación generados por Octave. No se compara aquí el preprocesado
alternativo skimage/farthest-point. La identidad de componentes se resolvió
mediante asignación global por SAD.

| Métrica Python frente a Octave | Resultado | Umbral 0,001 |
|---|---:|---|
| SAD máximo (rad), componente 7 | 0,001309545284 | Excede |
| RMSE global de endmembers | 0,000961700352 | Dentro |
| RMSE global de abundancias | 0,004907707362 | Excede |
| RMSE de abundancia máximo, componente 1 | 0,009769297173 | Excede |

Los RMSE de abundancias de **los nueve componentes** exceden 0,001:
0,00976930; 0,00263106; 0,00240086; 0,00235320; 0,00254897;
0,00631655; 0,00539586; 0,00413707; 0,00323734.
El resto de los SAD individuales queda por debajo del umbral. El informe
completo, incluidos ambos resultados frente a A/S conocidos, está en
[`equivalence.json`](../results/equivalence/equivalence.json).
**La equivalencia estricta a 0,001 no está demostrada: la prueba completa
falla y el script devuelve código 1.** No se ha aumentado el umbral para
ocultar el resultado.

Frente a la abundancia conocida del demo, RMSE Octave = 0,04416526 y
Python = 0,04391997. Esta proximidad de desempeño no demuestra equivalencia
numérica entre implementaciones ni validación bioóptica.

### Diagnóstico

Se volvió a ejecutar una iteración desde `octave_inputs.mat`, calculando
también `fcls(W0,seg.X_c)` antes de actualizar. La diferencia máxima en
esas soluciones FCLS de SciPy y Octave fue 0,08562299. El sistema aumentado
`[1e-5 W0; ones]` tiene condición 1,940935×10⁶. Tras una sola iteración,
RMSE de abundancias = 1,30684×10⁻⁵ y RMSE de endmembers = 2,20677×10⁻⁶
([diagnostic.json](../results/equivalence/diagnostic.json)). Esto identifica
una diferencia en NNLS que puede propagarse a los pesos espaciales; no
prueba que explique por sí sola toda la desviación a 100 iteraciones.
Se mantienen las ecuaciones publicadas y el NNLS de cada biblioteca.
Una réplica de las tolerancias de `lsqnonneg` sería trabajo adicional.

Octave emitió un aviso `svds` durante VCA: devolvió menos valores singulares
convergidos de los solicitados. No se oculta su efecto: las matrices que
produjo se compartieron con Python. No se distribuyen los MAT originales;
el harness permite regenerarlos desde el commit documentado.

## Escena costera y reproducción

PACE OCI, 27-abr-2024 17:45:13, productos OC_AOP y OC_BGC V3.2,
bbox [-74.5, 11.5, -71, 13.5], frente a La Guajira, Colombia.
Acceso NASA autenticado mediante earthaccess + HyperCoast. La búsqueda
acotada y los descartes por calidad quedan en `scene_selection.json`.

- 730 píxeles espectrales válidos; 719 pares válidos con chlor_a estándar.
- P=4, 25 superpíxeles solicitados, λ=0,3, semilla 17, 100 iteraciones.
- Parada por límite de iteraciones: no se afirma convergencia.
- Repetición desde `data/coastal_snapshot.npz`: diferencia absoluta máxima
  **0,0** en ambas matrices, comparada con la ejecución autenticada local.
- Spearman de fracciones contra log10(chlor_a), por componente:
  0,55094; −0,17819; −0,96677; 0,44276. No se calculan significancias
  suponiendo independencia espacial ni se asignan identidades a partir de ellas.
- Desviación máxima de la suma de abundancias respecto de uno: 0,00138460.

La figura `results/caribbean_replay/chlorophyll_sanity.png` y la variable
`chlor_a_reference` del NetCDF de repetición permiten inspeccionar la
comparación. Los factores son fracciones espectrales exploratorias, no
concentraciones ni identificaciones confirmadas de agua, fitoplancton,
CDOM o sedimentos. Se necesita una biblioteca espectral medida compatible
y validación independiente para ese uso.

El fixture público independiente de HyperCoast (golfo de México,
23-abr-2024, V1.0 NRT) pasó el test de lectura y ajuste con 124 píxeles
válidos en el recorte de humo. La ejecución de mapas usó 501 píxeles y
se detuvo por estancamiento a 52 iteraciones. SHA256 verificado en el loader.

## Límites de la entrega

Pruebas locales: 12 aprobadas (incluido el fixture PACE público); la
prueba Octave opcional, ejecutada por separado con cinco iteraciones,
**falló**: RMSE de abundancias = 0,00347641. Se conserva la aserción estricta
en pytest y CI para que el fallo siga visible. El informe completo usa
otra realización del preprocesado Octave; semillas iguales no garantizan
idéntica salida del eigensolver entre ejecuciones.

EMIT tiene un wrapper pero no una escena real validada. La imagen Docker
`sgsnmf-coastal:0.1.0` se construyó correctamente en Docker Desktop el
4-sep-2026; queda por registrar la ejecución del contenedor. Las pruebas
opcionales omitidas no se contabilizan como aprobadas. La prueba corta de
Octave en CI comprueba cinco iteraciones y no reemplaza el informe de 100.
El mismo script de entrada `scripts/reproduce.py` se ejecutó fuera de Docker
sin red ni intervención y generó los productos de los 730 píxeles y el
sanity-check de 719 pares a partir del snapshot incluido.

GitHub Actions ejecutó también las pruebas en Ubuntu: el trabajo `unit`
aprobó y `octave` falló. El test público de descarga no se ejecuta en cada
push y quedó omitido en ese run. Evidencia:
https://github.com/xxxxbotanicaxxxx-collab/sgsnmf-coastal/actions/runs/33897543213
