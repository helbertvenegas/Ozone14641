# Informe de avance: grupos de ozono para extensiones dobles de Ore (14641)

**Autores:** James Gómez y Helbert Venegas  
**Corte:** 22 de noviembre de 2026 · **Repositorio auditado:** versión 7

## Conclusión

La versión 7 permite **demostrar los tres casos puntuales Q(1), R(1) y V(1)** mediante aritmética exacta y una prueba de exhaustividad: una cota superior de automorfismos que fijan elementos centrales explícitos y una cota inferior proporcionada por la conjugación con elementos normales no nulos. En los tres casos el grupo es C4. **No se ha demostrado todavía la clasificación para todos los valores PI del twist h.**

## Resultados obtenidos y qué prueban

| Familia | h | Resultado | Certificado | Alcance |
|---|---|---|---|---|
| Q | 1 | Oz ≅ C4 | `qrv_ozone_h1.py` | Central bidegrado (1,2), 8 candidatos → 4, generador Ad(x1) |
| R | 1 | Oz ≅ C4 | `r1_normal_certificate.py` | Central bidegrado (1,1), 8 candidatos → 4, generador Ad(x2*y1+x1*y2) |
| V | 1 | Oz ≅ C4 | `v1_normal_certificate.py` | Cota superior 4, normal de bidegrado (1,2), cuatro potencias distintas |
| Q,R,V | h general, raíz de unidad | **Pendiente** | — | Falta teorema uniforme o casos separados por orden de h |

Todos estos certificados dependen de la corrección de las presentaciones algebraicas de entrada. Las relaciones de Q y V proceden de `Texto pegado (4).txt`, Tabla `tab:Mt`; las de R se obtuvieron de la presentación pública de Zhang–Zhang, https://4d-as.ncag.info/families/ore-r/. Antes de someter el artículo conviene cotejar las matrices con la fuente original ZZ09.

## Panorama de las 26 familias

- **14 monomiales:** B, E, J, K, L, M, N, O, P, S, T, U, W, Z. Hay teoremas generales del **orden** de Oz en el manuscrito de 18 de noviembre; el repositorio contiene auditorías y controles representativos, no certificados independientes exhaustivos para parámetros generales.
- **D, G:** G(p,f,1) cuenta con demostración general **condicionada** al teorema de generación del centro; el caso h general sigue sin cerrar. D(p,h) no está cerrado; en D(i,1) la reducción conocida deja 32 candidatos frente a los 16 reportados.
- **Q, R, V:** cerrados **solo para h=1** en la versión 7; restantes twists pendientes. Las filas con asterisco de la tabla no tienen rango central certificado.
- **C, F, I:** existen casos PI tabulados, pero el criterio PI general y el grupo de ozono uniforme siguen abiertos en los borradores.
- **A, H, X, Y:** excluidas del locus PI según el manuscrito; la razón y referencias deben figurar en el texto principal, no como resultados del programa.

## Qué falta para que el repositorio respalde una afirmación exhaustiva

1. Introducir las matrices Sigma de ZZ09 para las 26 familias **de manera independiente** de las presentaciones derivadas, y fijar convenciones de índices, signos y PBW.
2. Para Q, R, V, calcular Oz para h general, incluyendo h=-1, h=i y raíces de orden impar, y resolver las cuatro filas con asterisco en cuanto a centro/PI-grado.
3. Cerrar D para h general y G para h distinto de 1, comprobando exhaustividad en característica cero.
4. Para las 14 monomiales, comprobar la acción sobre el centro y los factores invariantes; distinguir teorema manuscrito de certificado de software.
5. Para C, F, I, completar primero el criterio PI general y después la clasificación de Oz.
6. Ejecutar el cuaderno íntegro en un entorno limpio (Colab y GitHub Actions); no confundir las salidas ejecutadas antiguas con resultados nuevos.

## Auditoría de integridad / límites

- El archivo `ozone_table.csv` recoge **73 filas** de ejemplos; solo se verificó su aritmética, no una prueba independiente de cada grupo.
- `certificate_report.json` recoge **280 comparaciones exactas** de dos implementaciones equivalentes de relaciones en las 14 familias monomiales: esos ensayos **no** prueban por sí mismos la exhaustividad de Oz.
- Las pruebas de Q, R y V para h=1 usan operaciones simbólicas exactas y no muestras aleatorias ni cuerpos finitos.
- El cuaderno principal conserva secciones históricas; la sección v7 llama a certificados nuevos. El archivo ejecutado antiguo **no debe presentarse** como versión ejecutada v7.
- **Apto para compartir como código suplementario de verificación parcial** siempre que el README y esta delimitación acompañen al paquete. **No apto para afirmar verificación exhaustiva de todas las familias PI.**

## Archivos que debe incluir el repositorio GitHub

`ozone_groups_verification.ipynb`, `ozone_table.csv`, `family_catalog.csv`, los módulos Python `general_quadratic.py`, `degree2_certificates.py`, `source_monomial.py`, `run_certificates.py`, `qrv_candidates.py`, `qrv_center_low.py`, `qrv_ozone_h1.py`, `r1_normal_certificate.py`, `v1_normal_certificate.py`, `dg_exact.py`, `qrv_dynamics.py`, `requirements.txt`, `test_certificates.py`, `test_qrv_h1.py`, `README.md`, `VERIFICATION_STATUS.md`, `INFORME_AVANCE.md`, `QRV_H1_COMPLETE.tex`, y los documentos de prueba de etapas anteriores. Se excluyen `__pycache__` y copias de kernel temporales.

## Ejecución rápida

```bash
python -m pip install -r requirements.txt
python qrv_ozone_h1.py
python r1_normal_certificate.py
python v1_normal_certificate.py
python -m unittest test_qrv_h1 -v
```
