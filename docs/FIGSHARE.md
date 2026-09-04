# GitHub + Figshare

Figshare es el archivo de publicación elegido por el usuario. No se creará
una cápsula CodeOcean del port. El enlace a la cápsula original de Wang se
mantiene exclusivamente como procedencia de su licencia MIT.

Preparar un depósito de tipo **Software**, con título:
«SGSNMF Coastal: Python port and reproducible coastal PACE workflow».
Los autores del depósito son los responsables del port; Wang et al.,
HyperCoast y earthaccess se incluyen como referencias, sin atribuirles
participación en este desarrollo. La identidad del autor se confirma con
el propietario antes de publicar.

Archivos: ZIP del código, wheel instalable, informe de equivalencia,
resultados costeros NetCDF/PNG/JSON y un snapshot pequeño del recorte
analizado cuando esté disponible. Los datos MATLAB originales y tokens
no se incluyen. Conservar LICENSE, THIRD_PARTY.md y el MIT original.
Seleccionar una licencia de software compatible con MIT; si la cuenta
solo ofrece licencias de datos, no aplicar CC0 a todo el port de forma
automática: mantener el depósito en borrador y resolver la licencia.

Referencias del depósito:
- URL real del repositorio GitHub y tag de la versión archivada.
- https://doi.org/10.1109/TGRS.2017.2724944
- https://doi.org/10.21105/joss.07025
- https://github.com/earthaccess-dev/earthaccess
- GranuleUR, colección y procesado NASA registrados en run.json.

Reservar un DOI en el borrador cuando esté disponible; añadirlo a README
y CITATION.cff antes de empaquetar la versión definitiva. La reserva no
equivale a publicación. Después de publicar, verificar la página pública,
descargas y enlace DOI. No afirmar que existe un DOI si aún no se asignó.

Reproducción local:

```sh
docker build -t sgsnmf-coastal .
docker run --rm -v "$PWD/results:/results" sgsnmf-coastal
```

También puede ejecutarse sin Docker, dentro del environment fijado:

```sh
RESULTS_DIR=results/reproduced SGSNMF_CACHE=.cache python scripts/reproduce.py
```

En PowerShell configurar las variables con `$env:`. Para NASA proporcionar
EARTHDATA_USERNAME y EARTHDATA_PASSWORD mediante el gestor de secretos del
entorno, además de SGSNMF_SCENE, SGSNMF_COLLECTION, SGSNMF_BBOX y opcionalmente
SGSNMF_CHL_SCENE. Nunca empaquetar esos valores. La reproducción pública
predeterminada no requiere credenciales NASA.

Guías oficiales:
- https://help.figshare.com/article/how-to-upload-and-publish-your-research
- https://info.figshare.com/user-guide/how-to-reserve-a-doi/
