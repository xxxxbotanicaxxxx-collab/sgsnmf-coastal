# SGSNMF Coastal

Port Python del optimizador de Wang et al. (2017), con acceso a PACE/EMIT
mediante HyperCoast + earthaccess, máscaras L2, CLI y salida reproducible
NetCDF/PNG/JSON. **Producto de investigación:** las fracciones espectrales
no son concentraciones de clorofila, CDOM o sedimentos.

La revisión MATLAB de diciembre de 2018 usa **gradiente proyectado**, no
actualización multiplicativa. `skimage.slic` es una alternativa de
preprocesamiento, no una réplica angular exacta. Consulte
[PORTING.md](docs/PORTING.md) para firmas, ecuaciones y diferencias.

**Validación numérica:** en el demo completo de 100 iteraciones, RMSE de
abundancias Python–Octave = 0,00490771 y SAD máximo = 0,00130955 rad.
Ambos exceden 0,001; esta versión no afirma equivalencia estricta.
[Informe y diagnóstico](docs/VALIDATION.md).

## Instalación

```sh
conda env create -f environment.yml
conda activate sgsnmf-coastal
```

Alternativa Python 3.12:

```sh
python -m venv .venv
# Activar .venv según el sistema operativo.
python -m pip install -r requirements.lock.txt
python -m pip install --no-deps -e .
```

El environment fija las dependencias principales; `requirements.lock.txt`
fija también las transitivas del entorno verificado, con marcadores para
paquetes exclusivos de Windows. El Dockerfile Linux está preparado; su
imagen `sgsnmf-coastal:0.1.0` se construyó correctamente con Docker Desktop
el 4-sep-2026. El contenedor completó el replay incluido y escribió los cinco
productos previstos en `/results`.

## Ejemplo completo sin credenciales ni descarga manual

La ejecución NASA verificada frente a La Guajira (27-abr-2024, V3.2) se
reproduce **sin red ni credenciales** con el recorte archivado:

```sh
sgsnmf-coastal run --snapshot data/coastal_snapshot.npz --p 4 \
  --segments 25 --max-iter 100 --seed 17 --output results/replayed
```

El recorte contiene 730 píxeles válidos y 719 pares con chlor_a estándar.
La reproducción local coincidió exactamente en endmembers y abundancias
con la ejecución autenticada. Para probar además la descarga automática
del fixture público independiente de HyperCoast:

```sh
sgsnmf-coastal run --public-demo --bbox -91 27 -88 29 --p 4 \
  --segments 25 --max-side 96 --max-iter 100 --seed 17 \
  --output results/public_pace
```

Descarga automática (~192 MB) del mismo fixture público usado por los tests
de HyperCoast: `PACE_OCI.20240423T184658.L2.OC_AOP.V1_0_0.NRT.nc`. Es un
producto histórico V1.0 de prueba; no representa el procesado NASA vigente.
El recorte es del golfo de México. El hash del archivo, versiones,
parámetros, escala, máscara y limitaciones quedan en `run.json`.

Salidas: `abundances.nc`, `abundance_maps.png`, `endmembers.png`, `run.json`
y `chlorophyll_sanity.png` cuando existe referencia válida.
Latitud/longitud son coordenadas 2D de la pasada; no se finge una rejilla
geográfica regular. Suma uno es una restricción suave, y su desviación se
reporta sin renormalizar las abundancias.

## Caribe colombiano con NASA

Crear Earthdata Login en https://urs.earthdata.nasa.gov/users/new y configurar
`EARTHDATA_USERNAME` y `EARTHDATA_PASSWORD` como secretos del entorno. No se
escriben credenciales en archivos ni se pasan por argumentos CLI. El código
usa `earthaccess.login(strategy='environment', persist=False)`.

```sh
sgsnmf-coastal search --bbox -77.5 7.5 -74 12 \
  --dates 2024-04-20 2024-04-30 --collection PACE_OCI_L2_AOP

sgsnmf-coastal run \
  --scene PACE_OCI.20240427T174513.L2.OC_AOP.V3_2.nc \
  --collection PACE_OCI_L2_AOP --bbox -74.5 11.5 -71 13.5 --p 4 \
  --segments 25 --seed 17 --max-iter 100 \
  --chlorophyll-scene PACE_OCI.20240427T174513.L2.OC_BGC.V3_2.nc \
  --output results/caribbean
```

El AOP y el BGC se leyeron mediante earthaccess + HyperCoast con
autenticación NASA; los resultados están en `results/caribbean`. Las
colecciones NRT no retienen necesariamente
granulados históricos: aquí se usa la colección reprocesada. El bbox de
búsqueda indica intersección, no garantiza cobertura válida de toda la costa.

`earthaccess.open` mantiene acceso HTTPS por rangos mientras HyperCoast
abre el archivo y se materializa el recorte. La capa remota de coordenadas
debe leerse para localizarlo; no se promete coste de red cero. La prueba
pública usa caché automática porque su mirror GitHub distribuye archivos.

Para EMIT: `--mission emit --collection EMITL2ARFL --scene <id_real>`;
aportar control de calidad acuático. El wrapper EMIT no tiene una prueba
de escena real en esta entrega. No llamar «agua» a todos los píxeles de
reflectancia terrestre; la ausencia de flags falla por defecto.

## Endmembers y referencia chlor_a

Se admite P=3 o 4. [DOMAIN.md](docs/DOMAIN.md) define los regímenes y el
formato de biblioteca medida para `--library`. Sin biblioteca, los factores
se numeran; no se inventan identificaciones bioópticas. `--chlorophyll-scene`
o `--chlorophyll-file` añade un OC_BGC estándar de la misma pasada y versión
para una comparación cualitativa. El ejemplo público AOP no incluye chlor_a:
su informe indica expresamente que ese sanity-check está pendiente.

## Equivalencia y tests

```sh
git clone https://github.com/Xinyu-Wang/SGSNMF_TGRS.git work/upstream
git -C work/upstream checkout e124bc76499bbbeebe10001b410a11b1d4945966
# GNU Octave y paquete image instalados, sin licencia MATLAB.
python scripts/equivalence.py --upstream work/upstream --octave octave-cli \
  --max-iter 100 --output results/equivalence
pytest -q
RUN_NETWORK_TESTS=1 pytest -q -m integration
SGSNMF_UPSTREAM=work/upstream OCTAVE=octave-cli pytest -q -m octave
```

En PowerShell usar `$env:RUN_NETWORK_TESTS='1'` y las otras variables de
forma equivalente. Los tests de red y Octave se omiten explícitamente si
faltan sus prerrequisitos; un skip no cuenta como validación. El harness
principal ejecuta el demo original completo de 100×100×221 con nueve
endmembers, comparte inicialización y segmentación entre lenguajes y
documenta toda diferencia >1e-3.
La prueba corta Octave de cinco iteraciones también falló localmente
(RMSE de abundancias 0,00347641); su aserción sigue activa en CI. Los
12 tests locales restantes, incluido el fixture PACE, aprobaron.

## GitHub, Figshare y reproducción

[Preparación de Figshare](docs/FIGSHARE.md). El script
`scripts/reproduce.py` genera `/results` sin interacción, con un fixture
archivado predeterminado y opción NASA por secretos. No depende de MATLAB.
GitHub contiene el código versionado. Figshare publicó el port con DOI
[10.6084/m9.figshare.33437689](https://doi.org/10.6084/m9.figshare.33437689)
y archivó esta versión y sus resultados bajo licencia MIT.

## Créditos, citas y licencia

- Wang, X., Zhong, Y., Zhang, L. y Xu, Y. (2017). *Spatial Group Sparsity
  Regularized Nonnegative Matrix Factorization for Hyperspectral Unmixing*.
  IEEE TGRS 55(11), 6287–6304. https://doi.org/10.1109/TGRS.2017.2724944
- Liu, B. y Wu, Q. (2024). *HyperCoast: A Python Package for Visualizing and
  Analyzing Hyperspectral Data in Coastal Environments*. JOSS 9(100), 7025.
  https://doi.org/10.21105/joss.07025
- earthaccess contributors, NASA Earthdata y comunidad:
  https://github.com/earthaccess-dev/earthaccess

MIT para el wrapper y port con conservación de los avisos originales.
El GitHub de Wang no tiene licencia explícita, pero la cápsula pública de
los autores contiene licencia MIT y los archivos auditados son idénticos.
[THIRD_PARTY.md](THIRD_PARTY.md) documenta commits, hashes y el aviso literal.
Se cita el algoritmo original; no se reclama su reautoría. Ver CITATION.cff.
