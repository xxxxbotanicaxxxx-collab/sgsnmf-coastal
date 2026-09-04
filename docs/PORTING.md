# Auditoría del port

Fuente: `Xinyu-Wang/SGSNMF_TGRS`, revisión
`e124bc76499bbbeebe10001b410a11b1d4945966`. Los archivos están realmente en
`code/Func/`, no directamente en `code/`. El README declara la revisión de
7-dic-2018, aunque los encabezados internos aún dicen 27-feb-2018.

## Contrato MATLAB real

```matlab
[W,H] = sgsnmf(para,seg)
seg = slic_HSI(im,k,m,seRadius,nItr)
```

`para.X`: L×N, bandas por píxeles; `para.W`: L×M, endmembers iniciales;
`para.H`: M×N, abundancias iniciales; `para.M`: rango M. Se requiere este
último campo aunque el encabezado no lo enumera. `para.Y` se usa en el demo,
pero no se lee en el optimizador. `para.lambda=0.3` lo establece el demo,
no `default_SGSNMF`. Los otros valores predeterminados son `tol=0.05`,
`maxiter=100`, `timelimit=600` segundos CPU, `verbose=1`, `print_iter=5`.
El retorno W tiene L×M y H, M×N. No se reordenan ni normalizan al terminar.

`im` es filas×columnas×bandas. `k` es el número deseado de superpíxeles;
`m` pondera la distancia espacial (el demo usa 0.5). `seRadius=1`,
`nItr=10` son opcionales. El encabezado enumera cinco salidas antiguas,
pero el código revisado devuelve **un struct seg**:

| Campo | Contenido |
|---|---|
| `X_c` | L×G, espectros medios de cada grupo |
| `P` | G, número efectivo de superpíxeles; no número de endmembers |
| `Cj` | N×1, confianza por píxel |
| `labels` | N×1, etiquetas desde 1, aplanamiento por columnas |
| `Sw` | separación S de la cuadrícula |

En Python `p` es el número de endmembers y `n_segments` el número solicitado
de superpíxeles. Las etiquetas comienzan en 0; -1 marca píxeles excluidos.
Todo el aplanamiento espacial usa `order='F'`.

## Ecuaciones y detalles conservados

La versión publicada **no usa actualizaciones multiplicativas**. Resuelve
subproblemas mediante gradiente proyectado, hasta 100 pasos internos y 20
intentos de búsqueda de paso, con α inicial 1 y β=0.1. El criterio es
`0.99 <grad,d> + 0.5 <(WᵀW)d,d> < 0`. La curvatura de la penalización no
se incluye en este criterio porque así está escrito el código original.

En cada iteración se calculan pesos de grupo
`w_p = 1 / (M² FCLS(W, X_c) + 1)`. FCLS resuelve NNLS con matrices
`[1e-5 W; 1]` y `[1e-5 X_c; 1]`, después limita a [eps,1]; no impone una
igualdad exacta de suma uno. La actualización de H usa matrices aumentadas
`[W;15]` y `[X;15]` y el gradiente adicional por píxel j
`lambda * Cj * w_p² * H_j / ||w_p * H_j||₂`.
La estructura espacial reside en pesos compartidos y confianza de grupo;
la norma en esta implementación se calcula sobre componentes de cada
píxel, no sobre todos los píxeles del grupo.

El historial de parada original comienza con `[SSE_inicial, 0]`; se conserva
internamente este cero artificial. Se detiene después de cinco iteraciones
consecutivas y veinte acumuladas sin una disminución SSE > 1e-4 (las cuatro
primeras reinician el contador consecutivo). La SSE es solamente el error
de reconstrucción, no el objetivo regularizado completo. El historial
exportado omite el cero artificial. El límite temporal es CPU, por lo que
se desactiva efectivamente para comparaciones deterministas.

## Diferencias deliberadas

1. **SLIC:** skimage usa distancia euclídea multibanda, otra inicialización y
   limpieza de conectividad. El original usa cuadrícula hexagonal y
   `sqrt(acos(coseno) + m² * distancia_espacial²/S²)`; no reinicia las
   distancias entre iteraciones. Su confianza usa `m`, no `m²`.
   `skimage.slic` no es un equivalente numérico directo. El adaptador conserva
   la fórmula de confianza y los promedios a partir de sus propias etiquetas;
   estima S con el tamaño medio de los grupos.
2. **Inicialización:** farthest-point determinista en Python en vez del VCA
   aleatorio original. No se confunde igualdad de semilla con igualdad de
   generadores aleatorios entre lenguajes.
3. **NNLS:** SciPy `nnls` frente a Octave `lsqnonneg`; el escalado 1e-5 puede
   hacer relevantes las tolerancias y el condicionamiento. Las diferencias
   se cuantifican en el informe, no se corrigen sustituyendo el original.
4. **Columnas nulas:** el cociente 0/0 del gradiente original se sustituye por
   un subgradiente cero. Se valida explícitamente este caso degenerado.
5. **Robustez:** entradas inválidas producen errores explicativos y máscaras
   excluyen datos ausentes. No se pintan abundancias de ceros sobre tierra.
6. **Métricas:** asignación global Hungarian por SAD, en vez del emparejamiento
   greedy del script original; la misma asignación se aplica a H. SAD se
   reporta en radianes y RMSE conserva la escala, sin ajustes post hoc.

## Demo y protocolo numérico

`data/F1_A_9.mat`: A de 221×9; `F1_S_9.mat`: S de 100×100×9. No contiene una
escena costera: son endmembers y mapas de abundancia usados para sintetizar
X=A·S. El demo agrega ruido gaussiano blanco de 20 dB y trunca a eps; solicita
`round(10000/36)=278` superpíxeles. VCA sobre medias de superpíxeles inicializa W,
y FCLS inicializa H.

`scripts/equivalence.py` ejecuta el original en GNU Octave con su paquete
image. Solo suprime el bloque que abre figuras de slic_HSI en una copia
temporal; conserva sus ecuaciones, VCA, NNLS y segmentación. Desactiva el
aviso repetido `Octave:possible-matlab-short-circuit-operator` para evitar
cientos de MB de mensajes dentro de los bucles; no cambia la operación `|`.
Exporta a MAT
el ruido realizado, X, W0, H0, seg y resultado. Python recibe exactamente
esas entradas. Esto aísla la equivalencia del **optimizador**, y no afirma
equivalencia de toda la cadena skimage/farthest-point. Se comparan W/H entre
lenguajes y ambos contra A/S, con umbral 1e-3 y salida no cero si se excede.

Los resultados efectivos y discrepancias se guardan en
`results/equivalence/equivalence.json`. El archivo MAT de referencia no se
incluye en Git. Volver a ejecutar el comando genera el artefacto completo.
