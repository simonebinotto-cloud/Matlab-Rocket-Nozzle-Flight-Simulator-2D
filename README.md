# 🚀 2D Rocket Flight & De Laval Nozzle MATLAB Simulator

A comprehensive numerical simulator developed in MATLAB for analyzing the flight dynamics of a launch vehicle / sounding rocket. The model couples 2D flight kinematics and dynamics (featuring a *Pitch-Over* maneuver and *Gravity Turn*) with a quasi-1D converging-diverging (de Laval) nozzle thermodynamic model, solving for normal shock wave position and operating regimes as ambient pressure changes with altitude.

---

## 📌 Key Features

- **De Laval Nozzle Model (Quasi-1D):** Dynamic calculation of engine performance ($P_2$, $V_2$, $\dot{m}$, Thrust $T$) as atmospheric pressure varies. Automatic handling of all nozzle regimes: subsonic, internal normal shock wave in the divergent section (solved via bisection), and supersonic/choked flow.
- **2D Flight Dynamics (3-DoF Point Mass):** Vectorial decomposition along horizontal and vertical axes ($x, h, v_x, v_h$).
- **Gravity Turn Maneuver:** Automated pitch control divided into 3 distinct phases (Pure vertical ascent, Forced pitch-over maneuver, and Vectorial gravity turn).
- **Standard Atmosphere (US Standard Atmosphere 1976):** Layered atmospheric model calculating temperature, static pressure, air density, and local speed of sound up to 47 km altitude.
- **Variable Aerodynamics:** Drag coefficient $C_d(\text{Mach})$ interpolated dynamically based on the instantaneous Mach number.
- **Dynamic Pressure & Max-Q Analysis:** Continuous monitoring of dynamic pressure $q = \frac{1}{2}\rho v^2$ with automatic identification of the maximum structural load point.

---

## 📐 Physical & Thermodynamic Model

### 1. Nozzle Thrust & Operating Regimes
Dynamic thrust is calculated at each integration step using the thrust equation:
$$T = \dot{m} V_2 + (P_2 - P_{\text{amb}}) A_2$$

When ambient pressure $P_{\text{amb}}$ falls within the over-expanded regime with an internal normal shock, the shock position $x_{\text{shock}}$ inside the divergent duct is solved numerically using Rankine-Hugoniot jump relations.

### 2. 2D Equations of Motion
$$\frac{dx}{dt} = v_x$$
$$\frac{dh}{dt} = v_h$$
$$\frac{dv_x}{dt} = \frac{T - D}{m(t)} \cos(\theta)$$
$$\frac{dv_h}{dt} = \frac{T - D}{m(t)} \sin(\theta) - g(h)$$

where $D = \frac{1}{2} \rho v^2 A C_d(\text{Mach})$ is the aerodynamic drag and $\theta$ represents the pitch angle.

---

## 📁 Project Structure

The repository contains the following progressive versions of the simulator:

- `nozzle_analysis.m`: Isolated script for thermodynamic analysis of the de Laval nozzle and regime transition under varying backpressure.
- `rocket_simulator_1D.m`: 1D atmospheric flight simulator (purely vertical trajectory).
- `rocket_simulator_2D.m`: **(Main Script)** Full 2D simulator featuring *Gravity Turn*, *Max-Q* analysis, and a comprehensive graphical dashboard.
