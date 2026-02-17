# -*- coding: utf-8 -*-
"""
Simulación gravitacional 2D: 4 cuerpos con animación
- Integrador: Verlet de velocidad (leapfrog)
- Suavizado gravitacional (epsilon) para evitar singularidades
- Centro de masa en reposo (marco del COM)

Sugerencias:
- Si notas inestabilidad numérica, reduce DT (p. ej. 0.005),
  aumenta SUBSTEPS (p. ej. 8) o incrementa EPS2 (p. ej. 1e-2).
- Cambia las condiciones iniciales en setup_four_bodies() para explorar otros regímenes.
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from collections import deque

# ====================== PARÁMETROS ======================
G = 1.0           # Constante gravitacional (unidades naturales)
DT = 0.01         # Paso de tiempo del integrador por sub-paso
SUBSTEPS = 5      # Sub-pasos por frame (mejora estabilidad visual)
EPS2 = 1e-3       # Suavizado gravitacional (epsilon^2)
TRAIL_LEN = 2400  # Longitud de la estela (historial)
WINDOW = 12.0     # Semialcance de la ventana de visualización (límites +/- WINDOW)
# ========================================================

class Body:
    """
    Cuerpo puntual con masa, posición, velocidad y color para graficar.
    """
    def __init__(self, x, y, vx, vy, m, color='tab:blue'):
        self.r = np.array([x, y], dtype=float)   # posición
        self.v = np.array([vx, vy], dtype=float) # velocidad
        self.m = float(m)                        # masa
        self.color = color
        self.trail = deque(maxlen=TRAIL_LEN)     # estela (historial)
        self.a = np.zeros(2, dtype=float)        # aceleración actual (cache)

def accelerations(bodies):
    """
    Calcula aceleraciones mutuas (gravedad Newtoniana con suavizado EPS2).
    Complejidad O(N^2), suficiente para N=4 en tiempo real.
    """
    # Reinicia aceleraciones
    for b in bodies:
        b.a[:] = 0.0

    N = len(bodies)
    for i in range(N - 1):
        ri = bodies[i].r
        mi = bodies[i].m
        for j in range(i + 1, N):
            rj = bodies[j].r
            mj = bodies[j].m
            dr = rj - ri
            dist2 = dr[0]*dr[0] + dr[1]*dr[1] + EPS2
            inv_r = 1.0 / math.sqrt(dist2)
            inv_r3 = inv_r / dist2
            # Fuerza gravitacional vectorial: F = G * mi * mj * dr * inv_r^3
            fvec = G * mi * mj * dr * inv_r3
            # a_i = F/mi ; a_j = -F/mj
            bodies[i].a += fvec / mi
            bodies[j].a -= fvec / mj

def step_verlet(bodies, dt):
    """
    Un paso de Verlet de velocidad (leapfrog):
      1) a(t)     = a[r(t)]
      2) v(t+dt/2)= v(t) + 0.5*dt*a(t)
      3) r(t+dt)  = r(t) + dt*v(t+dt/2)
      4) a(t+dt)  = a[r(t+dt)]
      5) v(t+dt)  = v(t+dt/2) + 0.5*dt*a(t+dt)
    """
    # a(t)
    accelerations(bodies)

    # v(t+dt/2) y r(t+dt)
    for b in bodies:
        v_half = b.v + 0.5 * dt * b.a
        b.r += dt * v_half
        b.v = v_half  # temporalmente guardamos v_half

    # a(t+dt)
    accelerations(bodies)

    # v(t+dt)
    for b in bodies:
        b.v += 0.5 * dt * b.a

def to_com_frame(bodies):
    """
    Traslada al marco del centro de masa: COM en (0,0) y velocidad del COM = 0.
    Mantiene el problema sin traslación ni deriva neta.
    """
    M = sum(b.m for b in bodies)
    r_com = sum(b.m * b.r for b in bodies) / M
    v_com = sum(b.m * b.v for b in bodies) / M
    for b in bodies:
        b.r -= r_com
        b.v -= v_com

def setup_four_bodies():
    """
    4 cuerpos: parte de una configuración de 3 cuerpos y añade uno más ligero.
    Ajusta posiciones/velocidades para explorar otras dinámicas.
    """
    # Masas
    m1, m2, m3, m4 = 1.0, 2.0, 3.0, 0.6

    # Posiciones iniciales (x, y)
    r1 = (-3.0,  2.4)
    r2 = ( 0.0,  0.0)
    r3 = ( 3.4,  2.8)
    r4 = ( 1.2, -3.2)

    # Velocidades iniciales (vx, vy)
    v1 = (-0.8,  0.0)
    v2 = ( 0.0,  0.0)
    v3 = ( 0.8,  0.0)
    v4 = ( 0.6,  0.35)

    b1 = Body(r1[0], r1[1], v1[0], v1[1], m1, color='tab:red')
    b2 = Body(r2[0], r2[1], v2[0], v2[1], m2, color='tab:gray')
    b3 = Body(r3[0], r3[1], v3[0], v3[1], m3, color='tab:blue')
    b4 = Body(r4[0], r4[1], v4[0], v4[1], m4, color='tab:green')
    bodies = [b1, b2, b3, b4]
    return bodies

# ======================= PREPARAR ESCENARIO ===========================
bodies = setup_four_bodies()
to_com_frame(bodies)

# ======================== CONFIGURAR GRÁFICA ==========================
fig, ax = plt.subplots(figsize=(8.0, 8.0))
ax.set_aspect('equal', 'box')
ax.set_xlim(-WINDOW, WINDOW)
ax.set_ylim(-WINDOW, WINDOW)
ax.grid(True, alpha=0.25)

# Discos y líneas para estelas
scatters = [ax.scatter(b.r[0], b.r[1], s=60, color=b.color, zorder=3) for b in bodies]
lines = [ax.plot([], [], color=b.color, lw=1.6, alpha=0.9)[0] for b in bodies]

def update(frame):
    # Integra varios sub-pasos por frame para mayor estabilidad visual
    for _ in range(SUBSTEPS):
        step_verlet(bodies, DT)

    # Actualiza estelas y marcadores
    for idx, b in enumerate(bodies):
        b.trail.append((b.r[0], b.r[1]))
        scatters[idx].set_offsets([b.r[0], b.r[1]])
        if len(b.trail) >= 2:
            xs, ys = zip(*b.trail)
            lines[idx].set_data(xs, ys)

    ax.set_title(f"Simulación 4 cuerpos | DT={DT}, substeps={SUBSTEPS}, suavizado={EPS2:.0e}")
    return scatters + lines

ani = FuncAnimation(fig, update, interval=16, blit=False)
plt.show()