# Aplicación del Método de Gradientes Conjugados Precondicionado (GCP) en Problemas Financieros

**Universidad Mayor de San Andrés (UMSA)**  
**Facultad de Ciencias Puras y Naturales — Carrera de Informática**  
**Materia:** Métodos Numéricos (INF-281)  
**Modalidad:** Individual  

---

## 1. Descripción General del Proyecto

Este repositorio contiene la implementación interactiva, rigurosa y autónoma del **Método de Gradientes Conjugados Precondicionado (GCP)** aplicado a un modelo de optimización de carteras y asignación presupuestaria interdepartamental ($n \ge 4$). 

El desarrollo se fundamenta en los principios algorítmicos expuestos por **Porfirio Suñagua S. (2020)** en su artículo científico *"Método de Gradientes Conjugados Precondicionado"* (Revista Boliviana de Matemática UMSA, Vol. 04, pp. 2–7). La solución contrasta el desempeño de **Cholesky (solución numérica de referencia directa)**, **Gradientes Conjugados simple (GC)** y **Gradientes Conjugados Precondicionado (GCP con Jacobi y SSOR)** a lo largo de tres escenarios operativos:
1. **Ideal** (Estable y bien condicionado).
2. **Bajo Estrés Financiero** (Perturbación por inflación y tasas).
3. **Mal Condicionado** (Severo mal escalamiento diagonal y alto acoplamiento).

---

## 2. Formulación Matemática del Modelo Financiero

### 2.1. Problema de Optimización Cuadrática
Una corporación tecnológica o de servicios financieros busca determinar el vector de asignación presupuestaria óptima $\mathbf{x} = [x_1, x_2, \dots, x_n]^\top \in \mathbb{R}^n$ (medido en millones de bolivianos / USD) entre $n=5$ unidades operativas.

La corporación busca minimizar una función cuadrática de costo de absorción, fricciones y riesgo cruzado, descontando los retornos proyectados:

$$\min_{\mathbf{x}} f(\mathbf{x}) = \frac{1}{2} \mathbf{x}^\top \mathbf{A} \mathbf{x} - \mathbf{b}^\top \mathbf{x}$$

Dado que la condición necesaria y suficiente de optimalidad de primer orden es $\nabla f(\mathbf{x}) = \mathbf{0}$:

$$\nabla f(\mathbf{x}) = \mathbf{A}\mathbf{x} - \mathbf{b} = \mathbf{0} \iff \mathbf{A}\mathbf{x} = \mathbf{b}$$

### 2.2. Significado Financiero de las Variables ($n=5$)
* **$x_1$ (I+D & Plataformas Cloud / IA):** Inversión en microservicios, analítica predictiva y modernización técnica.
* **$x_2$ (Infraestructura y Expansión Operativa):** Despliegue comercial, conectividad y redes logísticas.
* **$x_3$ (Ciberseguridad y Cumplimiento Regulatorio):** Mitigación de riesgos operacionales, auditorías y gobernanza de datos.
* **$x_4$ (Marketing Digital & Ventas):** Captación multicanal, fidelización de clientes y penetración de mercado.
* **$x_5$ (Fondo de Liquidez Estratégica & Tesorería):** Activos de reserva inmediata para cobertura frente a shocks externos.

### 2.3. Propiedades Requeridas de la Matriz $\mathbf{A}$
1. **Simetría ($\mathbf{A} = \mathbf{A}^\top$):** $a_{ij} = a_{ji}$ representa el impacto contractual, logístico y de fricción bilateral simétrica entre las unidades $i$ y $j$.
2. **Definida Positiva ($\mathbf{x}^\top \mathbf{A} \mathbf{x} > 0, \forall \mathbf{x} \neq \mathbf{0}$):** Garantiza que la función cuadrática $f(\mathbf{x})$ posea un mínimo global estricto y único. Se evalúa de manera automática mediante el **Criterio de Sylvester** (todos los menores principales estrictamente positivos: $\Delta_k = \det(\mathbf{A}_{1..k, 1..k}) > 0$).

---

## 3. Algoritmo de Suñagua (2020) e Implementación Numérica

### 3.1. Esquema Iterativo GCP
Dado el sistema $\mathbf{A}\mathbf{x} = \mathbf{b}$, vector inicial $\mathbf{x}_0 = \mathbf{0}$ y residuo $\mathbf{r}_0 = \mathbf{b} - \mathbf{A}\mathbf{x}_0$:

1. Resolver el sistema de precondicionamiento: $\mathbf{M} \mathbf{z}_0 = \mathbf{r}_0$.
2. Fijar dirección de descenso inicial: $\mathbf{p}_0 = \mathbf{z}_0$.
3. Para $k = 0, 1, 2, \dots, N_{\max}$ hasta que $\frac{\|\mathbf{r}_k\|_2}{\|\mathbf{b}\|_2} < \epsilon = 10^{-6}$:
   * $\alpha_k = \frac{\mathbf{r}_k^\top \mathbf{z}_k}{\mathbf{p}_k^\top \mathbf{A} \mathbf{p}_k}$
   * $\mathbf{x}_{k+1} = \mathbf{x}_k + \alpha_k \mathbf{p}_k$
   * $\mathbf{r}_{k+1} = \mathbf{r}_k - \alpha_k \mathbf{A} \mathbf{p}_k$
   * Evaluar convergencia: Si $\frac{\|\mathbf{r}_{k+1}\|_2}{\|\mathbf{b}\|_2} < 10^{-6}$, detener.
   * Resolver $\mathbf{M} \mathbf{z}_{k+1} = \mathbf{r}_{k+1}$
   * $\beta_k = \frac{\mathbf{r}_{k+1}^\top \mathbf{z}_{k+1}}{\mathbf{r}_k^\top \mathbf{z}_k}$
   * $\mathbf{p}_{k+1} = \mathbf{z}_{k+1} + \beta_k \mathbf{p}_k$

### 3.2. Precondicionadores Soportados
* **Jacobi (Diagonal):** $\mathbf{M} = \text{diag}(\mathbf{A}) \implies z_i = \frac{r_i}{a_{ii}}$. Recomendado formalmente para matrices con heterogeneidad de escalas en la diagonal principal.
* **SSOR (Symmetric Successive Over-Relaxation):** Para sistemas con fuerte correlación fuera de la diagonal.
* **Identidad ($\mathbf{M} = \mathbf{I}$):** Corresponde al algoritmo de Gradientes Conjugados clásico sin precondicionar, utilizado como contraste de rendimiento.

### 3.3. Método Directo de Referencia
Se utiliza la **Factorización de Cholesky** ($\mathbf{A} = \mathbf{L}\mathbf{L}^\top$) seguida de sustitución progresiva ($\mathbf{L}\mathbf{y} = \mathbf{b}$) y regresiva ($\mathbf{L}^\top \mathbf{x} = \mathbf{y}$), catalogada como solución numérica exacta de referencia.

---

## 4. Descripción de los Tres Escenarios Obligatorios

| Característica / Indicador | 1. Escenario Ideal | 2. Bajo Estrés Financiero | 3. Mal Condicionado |
| :--- | :--- | :--- | :--- |
| **Contexto Financiero** | Estabilidad macroeconómica; crecimiento orgánico equilibrado. | Shock inflacionario, encarecimiento de deuda y contracción de mercado. | Desbalance severo de escala entre proyectos macro y micro-iniciativas. |
| **Modificación Matemática** | Matriz $\mathbf{A}$ balanceada y diagonalmente dominante. | Incremento de costos ($a_{ii}$ y $a_{ij}$) y reorientación de $\mathbf{b}$ hacia liquidez. | Mal escalamiento diagonal masivo en $\mathbf{A}$ ($10^3$ vs $10^{-1}$). |
| **Número de Condición $\kappa(\mathbf{A})$** | $\kappa \approx 3.7$ | $\kappa \approx 8.6$ | $\kappa > 15,000$ |
| **Condición Precondicionada $\kappa(\mathbf{M}^{-1}\mathbf{A})$** | $\kappa \approx 2.4$ | $\kappa \approx 4.1$ | $\kappa \approx 35.8$ (Reducción dramática) |
| **Iteraciones GC (Simple)** | $\le 5$ iteraciones | $\le 7$ iteraciones | Decenas / Estancamiento numérico |
| **Iteraciones GCP (Jacobi)** | $\le 5$ iteraciones | $\le 6$ iteraciones | $\le 5$ iteraciones (Rápida convergencia) |

---

## 5. Arquitectura del Software

La aplicación interactiva está desarrollada como una **Single-Page Application (SPA)** autocontenida sin dependencias de servidor:
* **Estructura e Interfaz:** HTML5 semántico con [Tailwind CSS](https://tailwindcss.com/) para diseño responsivo moderno y accesible.
* **Visualización de Datos:** [Chart.js](https://www.chartjs.org/) para:
  * Gráfico semilogarítmico de convergencia residual: $\log_{10}(\|\mathbf{r}_k\|_2 / \|\mathbf{b}\|_2)$ vs iteración $k$.
  * Gráfico de barras de asignación óptima de capital ($\mathbf{x}_i$ en M BOB).
* **Motor Numérico:** JavaScript puro en el cliente (`index.html`), implementando:
  * Determinantes por eliminación Gaussiana para el Criterio de Sylvester.
  * Factorización de Cholesky.
  * Solvers GCP y GC.
  * Estimación espectral de números de condición ($\kappa(\mathbf{A})$ y $\kappa(\mathbf{M}^{-1}\mathbf{A})$).
  * Métricas de error componente por componente: absoluto ($E_{a,i}$), relativo ($E_{r,i}$) y porcentual ($E_{p,i}$).

---

## 6. Módulo de Pruebas de Calidad de Software (QA)

La aplicación integra una suite de pruebas automatizadas que verifica el comportamiento ante condiciones operativas normales y anómalas:

1. **CP-01 (Carga de Escenarios):** Ejecución de los escenarios preestablecidos verificando convergencia con residuo menor a $10^{-6}$ en $\le n$ pasos teóricos.
2. **CP-02 (Validación de Matriz No Válida):** Inserción deliberada de una matriz asimétrica o no definida positiva ($\Delta_k < 0$). El sistema captura la excepción, previene el bloqueo de la interfaz y alerta visualmente al usuario.
3. **CP-03 (Edición Dinámica):** Modificación en tiempo real de cualquier coeficiente de $\mathbf{A}$ o $\mathbf{b}$ con recálculo reactivo instantáneo.
4. **CP-04 (Aceleración por Precondicionador):** Validación de que las iteraciones de GCP con precondicionador de Jacobi sean estrictamente menores o iguales a las de GC simple en el escenario mal condicionado.

---

## 7. Instrucciones de Instalación y Ejecución Local

No se requiere ningún gestor de paquetes (Node, npm, pip) ni servidor web dedicado.

### 7.1. Clonar el Repositorio
```bash
git clone https://github.com/TU_USUARIO/gcp-finanzas-umsa.git
cd gcp-finanzas-umsa
```

### 7.2. Ejecutar la Aplicación
Basta con abrir el archivo `index.html` en cualquier navegador web moderno (Google Chrome, Firefox, Safari, Edge, Brave):
* **Linux:** `xdg-open index.html`
* **macOS:** `open index.html`
* **Windows:** Doble clic en `index.html` o `start index.html`

O bien, servirlo localmente mediante Python:
```bash
python3 -m http.server 8000
# Abrir en el navegador: http://localhost:8000
```

---

## 8. Referencias Bibliográficas

1. **Suñagua S., P. (2020).** *Método de Gradientes Conjugados Precondicionado*. Revista Boliviana de Matemática UMSA, Vol. 04, pp. 2–7.
2. **Golub, G. H., & Van Loan, C. F. (2013).** *Matrix Computations* (4th ed.). Johns Hopkins University Press.
3. **Saad, Y. (2003).** *Iterative Methods for Sparse Linear Systems* (2nd ed.). Society for Industrial and Applied Mathematics (SIAM).

---

## 9. Datos del Desafío y Enlaces

* **Estudiante:** Univ. Postulante Académico
* **Institución:** Universidad Mayor de San Andrés (UMSA) — Carrera de Informática
* **Despliegue Público (Web):** `https://TU_USUARIO.github.io/gcp-finanzas-umsa/`
* **Repositorio de Código:** `https://github.com/TU_USUARIO/gcp-finanzas-umsa`