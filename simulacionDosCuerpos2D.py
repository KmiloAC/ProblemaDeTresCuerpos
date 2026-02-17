# -*- coding: utf-8 -*-
"""
Simulación gravitacional 2D: SOLO 2 cuerpos con animación
- Integrador: Verlet de velocidad (leapfrog)
- Suavizado gravitacional (epsilon) para evitar singularidades
- Centro de masa en reposo (marco del COM)
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from collections import deque

# ====================== PARÁMETROS ======================
G = 1.0           # Constante gravitacional (unidades naturales)
DT = 0.01         # Paso de tiempo del integrador
SUBSTEPS = 5      # Sub-pasos por frame (mejora estabilidad visual)
EPS2 = 1e-3       # Suavizado gravitacional (epsilon^2)
TRAIL_LEN = 2000  # Longitud de la estela
WINDOW = 8.0      # Semialcance de la ventana de visualización (límites +/- WINDOW)

# Masas y separación inicial
m1, m2 = 1.0, 3.0
D = 6.0  # separación inicial total entre los dos cuerpos
# ========================================================

class Body:
    def __init__(self, x, y, vx, vy, m, color='tab:blue'):
        self.r = np.array([x, y], dtype=float)   # posición
        self.v = np.array([vx, vy], dtype=float) # velocidad
        self.m = float(m)
        self.color = color
        self.trail = deque(maxlen=TRAIL_LEN)
        self.a = np.zeros(2, dtype=float)

def accelerations(bodies):
    """
    Calcula aceleraciones mutuas (gravedad Newtoniana con suavizado EPS2).
    """
    for b in bodies:
        b.a[:] = 0.0

    # Solo 2 cuerpos, pero dejamos genérico
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
            fvec = G * mi * mj * dr * inv_r3
            bodies[i].a += fvec / mi
            bodies[j].a -= fvec / mj

def step_verlet(bodies, dt):
    """
    Un paso de Verlet de velocidad (leapfrog):
      v(t+dt/2) = v(t) + 0.5*dt*a(t)
      r(t+dt)   = r(t) + dt*v(t+dt/2)
      a(t+dt)   = a[r(t+dt)]
      v(t+dt)   = v(t+dt/2) + 0.5*dt*a(t+dt)
    """
    # a(t)
    accelerations(bodies)

    # v(t+dt/2) y r(t+dt)
    for b in bodies:
        v_half = b.v + 0.5 * dt * b.a
        b.r += dt * v_half
        b.v = v_half

    # a(t+dt)
    accelerations(bodies)

    # v(t+dt)
    for b in bodies:
        b.v += 0.5 * dt * b.a

def to_com_frame(bodies):
    """
    Traslada al marco del centro de masa: COM en (0,0) y velocidad COM = 0.
    """
    M = sum(b.m for b in bodies)
    r_com = sum(b.m * b.r for b in bodies) / M
    v_com = sum(b.m * b.v for b in bodies) / M
    for b in bodies:
        b.r -= r_com
        b.v -= v_com

def setup_two_bodies(m1=1.0, m2=3.0, D=6.0):
    """
    Construye dos cuerpos en órbita casi circular alrededor del baricentro.
    """
    M = m1 + m2
    # Posiciones sobre el eje x, simétricas respecto al COM (que quedará en 0)
    r1x = -D * (m2 / M)
    r2x =  D * (m1 / M)

    # Velocidades iniciales para órbita circular aproximada:
    # Órbita circular: ω = sqrt(G*M / D^3). Velocidad: v_i = ω * |r_i|
    omega = math.sqrt(G * M / (D**3))
    v1y =  omega * abs(r1x)   # sentido +y
    v2y = -omega * abs(r2x)   # sentido -y (opuesto)

    b1 = Body(r1x, 0.0, 0.0, v1y, m1, color='tab:orange')
    b2 = Body(r2x, 0.0, 0.0, v2y, m2, color='tab:blue')
    return [b1, b2]

# Preparar escenario de 2 cuerpos
bodies = setup_two_bodies(m1=m1, m2=m2, D=D)
to_com_frame(bodies)

# Animación
fig, ax = plt.subplots(figsize=(7, 7))
ax.set_aspect('equal', 'box')
ax.set_xlim(-WINDOW, WINDOW)
ax.set_ylim(-WINDOW, WINDOW)
ax.grid(True, alpha=0.25)

scatters = [ax.scatter(b.r[0], b.r[1], s=70, color=b.color, zorder=3) for b in bodies]
lines = [ax.plot([], [], color=b.color, lw=1.6, alpha=0.9)[0] for b in bodies]

def update(frame):
    for _ in range(SUBSTEPS):
        step_verlet(bodies, DT)

    for idx, b in enumerate(bodies):
        b.trail.append((b.r[0], b.r[1]))
        scatters[idx].set_offsets([b.r[0], b.r[1]])
        if len(b.trail) >= 2:
            xs, ys = zip(*b.trail)
            lines[idx].set_data(xs, ys)

    ax.set_title(f"Simulación 2 cuerpos | DT={DT}, substeps={SUBSTEPS}, suavizado={EPS2:.0e}")
    return scatters + lines

ani = FuncAnimation(fig, update, interval=16, blit=False)
plt.show()