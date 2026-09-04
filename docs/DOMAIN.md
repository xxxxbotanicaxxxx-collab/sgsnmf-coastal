# Interpretación bioóptica

Objetivo de dominio: tres o cuatro regímenes espectrales de agua clara,
agua enriquecida en fitoplancton, CDOM y, opcionalmente, sedimentos.
En Rrs, absorción y retrodispersión de estos constituyentes interactúan
de forma no lineal. Una descomposición lineal NMF no identifica de forma
única sus concentraciones ni demuestra que un componente sea «agua pura».
Fondos someros, aerosoles residuales y brillo solar también alteran espectros.

Por eso el modo ciego produce componentes numerados. Para investigar las
hipótesis se admite `--library espectros.csv`, con primera columna
`wavelength_nm` y tres o cuatro columnas de espectros medidos de agua clara,
fitoplancton, CDOM y sedimentos, en ese orden. Deben ser espectros de la
misma magnitud radiométrica y unidades (Rrs para PACE), no coeficientes de
absorción mezclados con reflectancia. Se interpola dentro del rango observado,
nunca se extrapola. La librería inicializa W, que sigue libre durante el
ajuste; al terminar se emparejan los factores por SAD y se conservan etiquetas
de **hipótesis**. No se incluye una biblioteca inventada ni se atribuyen
valores mg/m³ a las fracciones NMF.

PACE se usa con Rrs L2 visible 400–700 nm, banderas de calidad y una escala
global positiva. Espectros con algún valor negativo o ausente en las bandas
seleccionadas se excluyen completamente. Esto es conservador y puede sesgar
la selección hacia aguas más brillantes; `valid_mask`, número de píxeles y
banderas rechazadas quedan en resultados. No se suma un offset que cambie
la forma espectral. EMIT ofrece reflectancia terrestre; su uso acuático
requiere control atmosférico y máscara de agua externos. `read_emit` está
disponible, pero la validación pública de esta entrega se realiza con PACE.

El sanity-check acepta un producto estándar OC_BGC con `chlor_a`:
`--chlorophyll-file` o `--chlorophyll-scene`. Se exige la misma geolocalización
(tolerancia 1e-5 grados). Cada producto conserva su máscara de calidad:
un fallo BGC excluye ese valor de chlor_a del sanity-check, sin invalidar un
espectro AOP que pasó sus propios controles.
Se reporta Spearman por componente contra log10(chlor_a), con al menos veinte
pares, sin p-valores que ignoren dependencia espacial. No se elige ni calibra
el componente de fitoplancton usando el mismo producto que se compara.
Ausencia de producto BGC se reporta como `not_available`, nunca como aprobado.

Fuentes primarias:
- Wang et al. 2017: https://doi.org/10.1109/TGRS.2017.2724944
- NASA PACE BGC y chlor_a: https://pace.oceansciences.org/pace_bgc_oci.htm
- NASA bioóptica y limitaciones costeras: https://oceancolor.gsfc.nasa.gov/files/resources/docs/technical/atbd_mod19.pdf
