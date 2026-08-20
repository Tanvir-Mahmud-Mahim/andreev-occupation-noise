"""Tabulate the universal weak-coupling BCS gap function u(t) = Delta(T)/Delta0."""
import os
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad

A = 1.7639

def rhs(u, t):
    def integrand(x):
        E = np.sqrt(x * x + u * u)
        a = A * E / t
        return 1.0 / E / (np.exp(min(a, 500)) + 1.0)
    v, _ = quad(integrand, 0, 60, limit=300)
    return 2.0 * v

def solve_u(t):
    if t >= 0.9999:
        return 0.0
    if t < 0.08:
        return 1.0 - np.sqrt(2 * np.pi * t / A) * np.exp(-A / t)
    return brentq(lambda u: np.log(1.0 / u) - rhs(u, t), 1e-6,
                  1.0 - 1e-13, xtol=1e-15)

if __name__ == "__main__":
    ts = np.concatenate([np.linspace(0.002, 0.08, 40, endpoint=False),
                         np.linspace(0.08, 0.95, 180),
                         np.linspace(0.9505, 0.9995, 50)])
    us = np.array([solve_u(t) for t in ts])
    out = os.path.join(os.path.dirname(__file__), "..", "data")
    os.makedirs(out, exist_ok=True)
    np.savez(os.path.join(out, "bcs_gap_table.npz"), t=ts, u=us)
    print("BCS gap table written;",
          f"u(0.5)={solve_u(0.5):.6f} (expect 0.956887)")
