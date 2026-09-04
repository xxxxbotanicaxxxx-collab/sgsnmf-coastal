# Recorte de análisis PACE

`coastal_snapshot.npz` contiene únicamente el recorte preprocesado de la
pasada `PACE_OCI.20240427T174513.L2.OC_AOP.V3_2.nc` y la referencia colocada
`PACE_OCI.20240427T174513.L2.OC_BGC.V3_2.nc`. Fuente: NASA OB.DAAC,
colecciones PACE_OCI_L2_AOP y PACE_OCI_L2_BGC, procesado 3.2.

Región: [-74.5, 11.5, -71, 13.5] (oeste, sur, este, norte), frente a La Guajira.
Se obtuvo mediante streaming autenticado HTTPS de earthaccess y loaders
HyperCoast. No contiene credenciales, identificadores de usuario ni tokens.

El array `cube` contiene Rrs dividido por una escala global; `scale` restaura
sr^-1. Se conservan longitudes de onda, latitud/longitud 2D, máscara, chlor_a
L2 con su máscara propia y metadatos JSON. No usa pickle. Los píxeles
inválidos son NaN y no intervienen en la factorización. El hash de este
archivo se registra en el informe de replay.

La licencia del software no se impone a los productos NASA originales.
Este snapshot derivado se distribuye con atribución a NASA PACE/OB.DAAC;
conservar los identificadores de las fuentes y los metadatos al reutilizarlo.
Documentación de productos: https://pace.oceansciences.org/pace_bgc_oci.htm
