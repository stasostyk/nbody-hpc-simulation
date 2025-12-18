# Integrators

The integrator advances the system by updating particle **positions and velocities** based on gravitational accelerations.
Different integrators approximate these variables differently, affecting accuracy, stability, energy conservation, and computational cost.

---

### Euler
Updates position using the old velocity, then updates velocity using the current acceleration. It is the simplest to implement but unstable with strong energy drift.

$$
s^{n+1} = s^n + v^n \Delta t, \quad
v^{n+1} = v^n + a^n \Delta t
$$


---

### Symplectic Euler
Updates velocity first and then advances position using the new velocity. Same cost as Euler, but much better long-term stability.

$$
v^{n+1} = v^n + a^n \Delta t, \quad
s^{n+1} = s^n + v^{n+1} \Delta t
$$


It's good for cheap long simulations.

---

### Velocity Verlet
Splits the velocity update into two half-steps and recomputes acceleration after the position update. Second-order, symplectic, and energy-stable. And manages to have a very good balance between accuracy and stability

$$
v^{n+\frac12} \rightarrow s^{n+1} \rightarrow a^{n+1} \rightarrow v^{n+1}
$$


It's usually the default integrator** for N-body simulations.

---

### Runge–Kutta 4 (RK4)
Uses four slope evaluations per timestep to achieve high accuracy. Very accurate short-term but not symplectic and more expensive.

$$
k_1, k_2, k_3, k_4 \Rightarrow \text{4th-order update}
$$


It's best performance is for short runs and validation because of its high accuracy.


