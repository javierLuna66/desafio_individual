# Desafío: Modelar un problema real con ecuaciones diferenciales

* **Estudiante:** Javier Luna Rojas
* **Materia:** Métodos Numéricos
* **Fecha Límite:** 6 de noviembre de 2026

---

## 1. Descripción del Problema y del Modelo Matemático

### Contexto Comunitario
Este proyecto aborda la crisis de sobrepoblación en los refugios temporales de rescate animal ubicados en el Distrito Municipal 3 de la ciudad de **El Alto**. El objetivo del modelo es ofrecer una herramienta cuantitativa para los gestores y voluntarios comunitarios, permitiendo proyectar en cuántas semanas se logrará reducir el hacinamiento inicial desde un estado crítico hasta una capacidad operativa sostenible, considerando el flujo constante de nuevos rescates y la efectividad de las campañas de adopción.

### Formulación de la EDO
El balance dinámico de la población canina $P$ en función del tiempo $t$ se rige por la ley de conservación: la razón de cambio poblacional es igual a la tasa de ingreso por rescates menos la tasa de salida por adopciones y reencuentros familiares.

$$\frac{dP}{dt} = I - kP, \quad P(0) = P_0$$

#### Variables y Parámetros del Sistema
| Símbolo | Descripción | Unidad | Valor Base |
| :---: | :--- | :---: | :---: |
| $t$ | Variable independiente: Tiempo transcurrido de campaña | $\text{semanas}$ | $t \in [0, 10]$ |
| $P(t)$ | Variable dependiente: Población de perros en el refugio | $\text{perros}$ | $P(0) = 80$ |
| $I$ | Tasa constante de ingreso por nuevos rescates | $\text{perros/semana}$ | $6.0$ |
| $k$ | Tasa proporcional de adopción y reencuentro | $\text{semana}^{-1}$ | $0.20$ ($20\%$) |
| $P_{\text{meta}}$ | Umbral de capacidad máxima sostenible del albergue | $\text{perros}$ | $45$ |

### Solución Analítica Exacta
Resolviendo la EDO lineal de primer orden mediante el factor integrante $\mu(t) = e^{kt}$, obtenemos la solución analítica exacta:

$$P(t) = \frac{I}{k} + \left(P_0 - \frac{I}{k}\right)e^{-kt}$$

Sustituyendo los parámetros base ($P_0 = 80$, $I = 6$, $k = 0.20$):

$$P(t) = 30 + 50e^{-0.20t}$$

* **Equilibrio asintótico ($t \to \infty$):** $P_{\text{eq}} = \frac{I}{k} = 30\text{ perros}$.
* **Tiempo exacto para alcanzar la meta ($P(t) = 45$):**
  $$t_{\text{meta}} = \frac{\ln\left(\frac{45 - 30}{80 - 30}\right)}{-0.20} = \frac{\ln(0.30)}{-0.20} \approx 6.0199\text{ semanas}$$

---

## 2. Métodos Numéricos Implementados

Definiendo $f(t, P) = I - kP$ y un tamaño de paso $h$ tal que $t_{i+1} = t_i + h$:

1. **Método de Heun (Runge-Kutta de 2.° Orden - RK2):**
   * Predictor: $k_1 = f(t_i, P_i)$
   * Corrector: $k_2 = f(t_i + h, P_i + h \cdot k_1)$
   * Actualización: $P_{i+1} = P_i + \frac{h}{2}(k_1 + k_2)$
   * Orden de error global: $O(h^2)$

2. **Método Clásico de Runge-Kutta de 4.° Orden (RK4):**
   * $k_1 = f(t_i, P_i)$
   * $k_2 = f\left(t_i + \frac{h}{2}, P_i + \frac{h}{2}k_1\right)$
   * $k_3 = f\left(t_i + \frac{h}{2}, P_i + \frac{h}{2}k_2\right)$
   * $k_4 = f(t_i + h, P_i + h \cdot k_3)$
   * Actualización: $P_{i+1} = P_i + \frac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4)$
   * Orden de error global: $O(h^4)$

3. **Métricas de Error:**
   * **Error Absoluto:** $E_a = |P_{\text{exacto}}(t_i) - P_{\text{numérico}, i}|$
   * **Error Relativo Porcentual:** $\varepsilon_r(\%) = \left|\frac{P_{\text{exacto}}(t_i) - P_{\text{numérico}, i}}{P_{\text{exacto}}(t_i)}\right| \times 100\%$

---

## 3. Arquitectura del Sistema Web

La aplicación web ha sido desarrollada bajo una arquitectura **100% Client-Side** en un archivo autocontenido (`index.html`), garantizando portabilidad y ejecución inmediata sin dependencias de backend:
* **HTML5 y CSS3:** Diseño responsivo con paleta de colores de alto contraste y grillas adaptables.
* **JavaScript (ES6):** Motor de cálculo iterativo, validación robusta de entradas y gestión de eventos en tiempo real.
* **Chart.js (v4):** Librería externa cargada por CDN para la representación gráfica interactiva de las trayectorias poblacionales y de los errores relativos.

---

## 4. Instrucciones de Ejecución y Despliegue

### Opción A: Ejecución en Línea (GitHub Pages)
1. Accede a la URL pública generada en el despliegue del repositorio:
   * **URL de la Aplicación:** `https://<tu-usuario>.github.io/desafio-edo-hogar-canino/`

### Opción B: Ejecución Local
1. Clona este repositorio en tu equipo:
   ```bash
   git clone https://github.com/tu-usuario/desafio-edo-hogar-canino.git
   cd desafio-edo-hogar-canino
   ```
2. Abre el archivo `index.html` directamente haciendo doble clic sobre él en cualquier navegador web moderno (Google Chrome, Mozilla Firefox, Microsoft Edge, Brave).

---
*Desafío individual desarrollado por **Javier Luna Rojas** para la materia de Métodos Numéricos.*