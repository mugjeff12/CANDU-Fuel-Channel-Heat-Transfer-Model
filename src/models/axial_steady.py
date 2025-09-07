# src/models/axial_steady.py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from src.props.props_water import water_props_T  # T-dependent (educational fits)

# ---------------- Base inputs ----------------
L = 6.0        # m  (channel length)
D_h = 0.012    # m  (hydraulic diameter of equivalent subchannel)
T_in = 560.0   # K  (inlet coolant temperature)
P = 10e6       # Pa (channel pressure)
m_dot = 2.0    # kg/s (mass flow rate)

# Choose heat rate (use 60e3 for "cooler" demo; 350e3 for "hot" base)
q_lin = 25e3  # W/m

# Pipe roughness for Churchill friction (m): ~5e-6 smooth metal
eps = 0

# Discretization
N = 200
x = np.linspace(0, L, N + 1)
dx = x[1] - x[0]

# Geometry
P_wet = np.pi * D_h
A_flow = np.pi * (D_h / 2) ** 2
q_pp = q_lin / P_wet  # W/m^2

# ---------------- Arrays ----------------
T_c = np.zeros_like(x); T_c[0] = T_in
T_w = np.zeros_like(x)
h   = np.zeros_like(x)
Re  = np.zeros_like(x)
Nu  = np.zeros_like(x)
f   = np.zeros_like(x)
dP  = np.zeros_like(x)  # cumulative Pa

# ---------------- Helpers ----------------
def h_from_props(rho, mu, k, cp, D_h, m_dot, A_flow):
    """Dittus–Boelter heat transfer using local props."""
    v  = m_dot / (rho * A_flow)
    Re = rho * v * D_h / mu
    Pr = cp * mu / k
    Nu = 0.023 * (Re**0.8) * (Pr**0.4)
    h  = Nu * k / D_h
    return h, Re, Nu, v

def f_churchill(Re, D_h, eps):
    """
    Churchill (1977) friction factor, valid for all Re and roughness.
    eps: absolute roughness [m]
    """
    if Re < 1e-12:
        return 0.0
    A = (2.457 * np.log((7.0/Re)**0.9 + 0.27*(eps/D_h)))**16
    B = (37530.0 / Re)**16
    return 8.0 * ((8.0/Re)**12 + 1.0/((A + B)**1.5))**(1.0/12.0)

def tsat_water(P_pa: float) -> float:
    """
    Saturation temperature for water at pressure P (Pa).
    Tries CoolProp; else returns an educational placeholder near 10 MPa.
    """
    try:
        import CoolProp.CoolProp as CP
        return CP.PropsSI("T", "P", P_pa, "Q", 0, "Water")
    except Exception:
        return 584.0  # ~Tsat at 10 MPa (K), for visualization only

# ---------------- Axial march with film-T iteration ----------------
# Inlet film properties & initial wall estimate
props0 = water_props_T(T_in, P)
h_guess, Re0, Nu0, v0 = h_from_props(
    rho=props0["rho"], mu=props0["mu"], k=props0["k"], cp=props0["cp"],
    D_h=D_h, m_dot=m_dot, A_flow=A_flow
)
T_w[0] = T_c[0] + q_pp / h_guess
h[0], Re[0], Nu[0] = h_guess, Re0, Nu0
f[0] = f_churchill(Re[0], D_h, eps)

for i in range(N):
    # 1) Bulk energy balance (cp at current film estimate)
    props_i = water_props_T(0.5*(T_c[i] + T_w[i]), P)
    T_c[i+1] = T_c[i] + (q_lin / (m_dot * props_i["cp"])) * dx

    # 2) Film-temperature fixed-point iteration at node i+1
    T_c_loc = T_c[i+1]
    T_w_loc = T_c_loc + q_pp / max(h[i], 1e-3)  # start from previous h
    pr = props_i
    for _ in range(3):  # few passes are sufficient here
        T_f = 0.5 * (T_c_loc + T_w_loc)
        pr  = water_props_T(T_f, P)
        h_loc, Re_loc, Nu_loc, v_loc = h_from_props(
            rho=pr["rho"], mu=pr["mu"], k=pr["k"], cp=pr["cp"],
            D_h=D_h, m_dot=m_dot, A_flow=A_flow
        )
        T_w_loc = T_c_loc + q_pp / h_loc

    # Store locals
    T_w[i+1] = T_w_loc
    h[i+1]   = h_loc
    Re[i+1]  = Re_loc
    Nu[i+1]  = Nu_loc
    f[i+1]   = f_churchill(Re_loc, D_h, eps)

    # 3) Local pressure gradient with updated rho, v
    dPdx = f[i+1] * (pr["rho"] * v_loc**2 / 2.0) / D_h
    dP[i+1] = dP[i] + dPdx * dx

# ---------------- Margin to boiling ----------------
T_sat = tsat_water(P)                 # K (assume constant along channel)
DeltaT_sat = T_sat - T_w              # K
min_margin = float(DeltaT_sat.min())
ix_min = int(np.argmin(DeltaT_sat))
x_min = float(x[ix_min])

# ---------------- Output paths ----------------
base = Path(__file__).resolve().parents[2]
data_dir = base / "data"; data_dir.mkdir(exist_ok=True, parents=True)
fig_dir  = base / "figures"; fig_dir.mkdir(exist_ok=True, parents=True)

# ---------------- Save CSV ----------------
out = pd.DataFrame({
    "x_m": x,
    "T_c_K": T_c,
    "T_w_K": T_w,
    "h_W_m2K": h,
    "Re": Re,
    "Nu": Nu,
    "f": f,
    "dP_cum_Pa": dP,
    "T_sat_K": np.full_like(x, T_sat, dtype=float),
    "DeltaT_sat_K": DeltaT_sat
})
out.to_csv(data_dir / "axial_steady_base.csv", index=False)

# ---------------- Plots ----------------
# 1) Temperatures + Tsat
plt.figure()
plt.plot(x, T_c, label="Coolant T_c", linewidth=2)
plt.plot(x, T_w, label="Wall T_w", linewidth=2, linestyle="--")
plt.plot(x, np.full_like(x, T_sat), label="T_sat (≈)", linewidth=1.5, linestyle="-.")
plt.xlabel("Axial position x (m)")
plt.ylabel("Temperature (K)")
plt.title("Coolant, Wall, and Saturation Temperature vs x")
plt.grid(True, linestyle=":")
plt.legend()
plt.savefig(fig_dir / "Tx_coolant_wall_Tsat.png", dpi=200, bbox_inches="tight")

# 2) Cumulative pressure drop
plt.figure()
plt.plot(x, dP/1e6, linewidth=2)
plt.xlabel("Axial position x (m)")
plt.ylabel("ΔP (MPa)")
plt.title("Cumulative Pressure Drop vs x (Churchill friction)")
plt.grid(True, linestyle=":")
plt.savefig(fig_dir / "dP_x_filmT.png", dpi=200, bbox_inches="tight")

# 3) Margin to boiling with min annotation
plt.figure()
plt.plot(x, DeltaT_sat, linewidth=2, label="ΔT_sat = T_sat - T_w")
plt.axhline(0, linestyle="--", linewidth=1)
plt.scatter([x_min], [min_margin], zorder=5)
plt.annotate(f"min {min_margin:.1f} K at x={x_min:.2f} m",
             (x_min, min_margin),
             textcoords="offset points", xytext=(8,8))
plt.xlabel("Axial position x (m)")
plt.ylabel("ΔT_sat (K)")
plt.title("Margin to Boiling vs x")
plt.grid(True, linestyle=":")
plt.legend()
plt.savefig(fig_dir / "margin_to_boiling.png", dpi=200, bbox_inches="tight")

# ---------------- Console summary ----------------
print("=== Axial Steady with Film-T + Churchill + Boiling Margin ===")
print(f"T_out = {T_c[-1]:.2f} K, T_w,out = {T_w[-1]:.2f} K, T_sat ≈ {T_sat:.2f} K")
print(f"h_in = {h[0]:.0f} W/m^2-K, h_out = {h[-1]:.0f} W/m^2-K")
print(f"Re range ~ {Re[0]:.2e} → {Re[-1]:.2e}")
print(f"Friction: Churchill, eps = {eps:g} m")
print(f"Total ΔP = {dP[-1]/1e6:.3f} MPa over L = {L} m")
print(f"Min ΔT_sat = {min_margin:.1f} K at x = {x_min:.2f} m")
print(f"CSV:    {data_dir / 'axial_steady_base.csv'}")
print(f"Figures:{fig_dir / 'Tx_coolant_wall_Tsat.png'}, {fig_dir / 'dP_x_filmT.png'}, {fig_dir / 'margin_to_boiling.png'}")

