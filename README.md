# Simulación del Problema de N Cuerpos (2, 3 y 4 cuerpos)
## Medición de Complejidad con Entropía de Shannon

Realizado por: Camilo Andres Lopez Contreras y Juan Jose Marquez Villareal

---

# 1. Objetivo

Este proyecto implementa simulaciones 2D y 3D de sistemas gravitacionales de **2, 3 y 4 cuerpos**, utilizando integración numérica estable y evaluando la **complejidad del sistema** mediante **Entropía de Shannon**

**“Para cada situación (2, 3 y 4 cuerpos), medimos la complejidad del sistema revisando el concepto de Entropía de Shannon.”**

---

# 2. Fundamentos Físicos

## 2.1 Ley de Gravitación Universal de Newton

La fuerza entre dos masas puntuales \(m_1\) y \(m_2\), separadas por una distancia \(r\), está dada por:

\[
F = G \frac{m_1 m_2}{r^2}
\]

En la simulación se utiliza \(G = 1.0\).

Las componentes en 2D son:

\[
F_x = F \frac{\Delta x}{r}, \qquad
F_y = F \frac{\Delta y}{r}, \qquad
r = \sqrt{(\Delta x)^2 + (\Delta y)^2}
\]

donde \(\Delta x = x_2 - x_1\) y \(\Delta y = y_2 - y_1\).

---

## 2.2 Segunda Ley de Newton

La aceleración de cada cuerpo se obtiene mediante:

\[
\vec{a} = \frac{\vec{F}_\text{neta}}{m}
\]

Cada cuerpo experimenta la suma de todas las fuerzas ejercidas por los demás cuerpos del sistema.

---

## 2.3 Suavizado Gravitacional

Para evitar singularidades numéricas cuando \(r \to 0\), se usa suavizado:

\[
r^2 \rightarrow r^2 + \varepsilon^2
\]

donde \(\varepsilon^2\) es pequeño.

---

# 3. Integración Numérica y Marco de Referencia

## 3.1 Integración con Verlet de Velocidad (Leapfrog)

Se utiliza el integrador **Verlet de velocidad**, que presenta mejor estabilidad y conservación de energía que Euler.

Esquema:

\[
\vec{v}\!\left(t+\tfrac{\Delta t}{2}\right)
= \vec{v}(t) + \tfrac{1}{2}\Delta t\,\vec{a}(t)
\]

\[
\vec{r}(t+\Delta t)
= \vec{r}(t) + \Delta t\,\vec{v}\!\left(t+\tfrac{\Delta t}{2}\right)
\]

Recalcular aceleraciones:

\[
\vec{v}(t+\Delta t)
= \vec{v}\!\left(t+\tfrac{\Delta t}{2}\right)
+ \tfrac{1}{2}\Delta t\,\vec{a}(t+\Delta t)
\]

---

## 3.2 Marco del Centro de Masa (COM)

Para evitar traslaciones globales del sistema, se lo transforma al marco donde:

- El COM está en el origen  
- Su velocidad es cero  

Esto mejora estabilidad y visualización.

---

# 4. Medición de Complejidad: Entropía de Shannon

## 4.1 Definición

Sea \(x(t)\) una serie temporal escalar. A partir de un histograma con probabilidades \(p_i\), la Entropía de Shannon se define como:

\[
H = -\sum_{i=1}^{N} p_i \log_2(p_i)
\]

La versión normalizada en \([0,1]\):

\[
H_{\text{norm}} = \frac{H}{\log_2(N)}
\]

---

## 4.2 Señal Utilizada

Para comparar entre 2, 3 y 4 cuerpos se usa la misma serie:

\[
d(t) = \left\lVert \vec{r}_0(t) - \vec{r}_{\text{COM}}(t) \right\rVert
\]

Es decir, la **distancia del cuerpo 0 al COM**.  
Sobre esta señal se calcula la entropía normalizada.

---

## 4.3 Interpretación

- \(H_{\text{norm}} \approx 0\): poca variabilidad en la señal.  
- \(H_{\text{norm}} \approx 1\): alta variabilidad y dispersión.  

**Nota:**  
Shannon mide **dispersión de valores**, no periodicidad ni caos dinámico.  
Un sistema periódico puede arrojar valores altos si la señal cubre un rango amplio.

# 6. Requisitos y Ejecución

## 6.1 Dependencias

```bash
pip install requirements.txt
python "Nombre-del-archivo.py"