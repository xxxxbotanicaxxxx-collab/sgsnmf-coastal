# Procedencia y licencias

SGSNMF no es un algoritmo nuevo de este proyecto. El núcleo Python adapta
el trabajo de Xinyu Wang, Yanfei Zhong, Liangpei Zhang y Yanyan Xu (2017).
Los encabezados MATLAB atribuyen Copyright (C) 2018 a Xinyu Wang y Yanfei
Zhong, Wuhan University, «All rights reserved».

El repositorio https://github.com/Xinyu-Wang/SGSNMF_TGRS, commit
`e124bc76499bbbeebe10001b410a11b1d4945966`, no contiene una licencia explícita.
Sin embargo, la cápsula pública de los autores
https://codeocean.com/capsule/2018026/tree/v2 sí publica `/code/LICENSE` MIT.
Se clonó desde la URL expuesta por su interfaz:
`https://git.codeocean.com/capsule-2018026.git`, commit
`1c6500939e6b7da2e5290a28b3936e0a42628cda`.

Se comprobó identidad SHA-256 de los dos archivos entre ambos checkouts:

| Archivo | SHA-256 de los checkouts locales |
|---|---|
| code/Func/sgsnmf.m | cd1241e35e22b1afa74ca9eff4474da0104547d7f99dca50db9a7b8fd4049ff1 |
| code/Func/slic_HSI.m | b70eda40d1899214f6807cae885b74aedc5665a237f23841aca98201bb0de78b |

Los hashes anteriores incluyen los finales de línea del checkout Windows;
la igualdad fue evaluada en la misma máquina. Se conserva el archivo MIT
literal en `licenses/SGSNMF-CodeOcean-MIT.txt`, incluidos sus marcadores
`[year] [fullname]`. No se rellenan ni se ocultan los marcadores originales.
La autoría original se mantiene explícita en este archivo y en LICENSE.

La licencia MIT del wrapper no cambia los derechos sobre datos externos.
Los MAT originales se obtienen del repositorio para la validación y no se
redistribuyen con el paquete. Los archivos de referencia Octave, que incluyen
esos datos, se excluyen de Git. La escena PACE pública procede del fixture de
HyperCoast y se almacena automáticamente en caché; tampoco se redistribuye.

HyperCoast: Bingqing Liu y Qiusheng Wu, JOSS 2024, MIT (algunas partes EMIT
conservan además avisos Apache-2.0). earthaccess: comunidad earthaccess,
NASA Earthdata y colaboradores; dependencia externa con sus propios avisos.
No se copian parsers HDF5/NetCDF ni credenciales.
