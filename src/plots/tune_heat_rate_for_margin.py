import numpy as np
from pathlib import Path
from src.props.props_water import water_props_T
import matplotlib.pyplot as plt

# --- copy of model core (functionized) ---
def h_from_props(rho, mu, k, cp, D_h, m_dot, A_flow):
    v  = m_dot / (rho * A_flow)
    Re = rho * v * D_h / mu
    Pr = cp * mu / k
    Nu = 0.023 * (Re**0.8) * (Pr**0.4)
    h  = Nu * k / D_h
    return h, Re, Nu, v

def f_churchill(Re, D_h, eps):
    if Re < 1e-12: return 0.0
    A = (2.457 * np.log((7.0/Re)**0.9 + 0.27*(eps/D_h)))**16
    B = (37530.0 / Re)**16
    return 8.0 * ((8.0/Re)**12 + 1.0/((A + B)**1.5))**(1.0/12.0)

def tsat_water(P_pa):
    try:
        import CoolProp.CoolProp as CP
        return CP.PropsSI("T", "P", P_pa, "Q", 0, "Water")
    except Exception:
        return 584.0

def run_channel(L, D_h, T_in, P, m_dot, q_lin, eps=0.0, N=200):
    x = np.linspace(0, L, N+1); dx = x[1]-x[0]
    P_wet = np.pi * D_h
    A_flow = np.pi * (D_h/2)**2
    q_pp = q_lin / P_wet

    T_c = np.zeros_like(x); T_c[0] = T_in
    T_w = np.zeros_like(x)
    h   = np.zeros_like(x)
    dP  = np.zeros_like(x)

    props0 = water_props_T(T_in, P)
    h0, Re0, Nu0, v0 = h_from_props(**{k:props0[k] for k in ["rho","mu","k","cp"]},
                                    D_h=D_h, m_dot=m_dot, A_flow=A_flow)
    T_w[0] = T_c[0] + q_pp / h0
    h[0]   = h0

    for i in range(N):
        props_i = water_props_T(0.5*(T_c[i]+T_w[i]), P)
        T_c[i+1] = T_c[i] + (q_lin/(m_dot*props_i["cp"])) * dx

        T_c_loc = T_c[i+1]
        T_w_loc = T_c_loc + q_pp / max(h[i], 1e-3)
        pr = props_i
        for _ in range(3):
            T_f = 0.5*(T_c_loc + T_w_loc)
            pr  = water_props_T(T_f, P)
            h_loc, Re_loc, Nu_loc, v_loc = h_from_props(**{k:pr[k] for k in ["rho","mu","k","cp"]},
                                                        D_h=D_h, m_dot=m_dot, A_flow=A_flow)
            T_w_loc = T_c_loc + q_pp / h_loc

        T_w[i+1] = T_w_loc
        h[i+1]   = h_loc
        f = f_churchill(Re_loc, D_h, eps)
        dPdx = f * (pr["rho"] * v_loc**2 / 2.0) / D_h
        dP[i+1] = dP[i] + dPdx * dx

    T_sat = tsat_water(P)
    margin = T_sat - T_w
    return dict(x=x, T_c=T_c, T_w=T_w, dP=dP, margin=margin, T_sat=T_sat)

# --- user targets/base geometry ---
L=6.0; D_h=0.012; T_in=560.0; P=10e6; m_dot=2.0; eps=0.0
target_margin = 10.0  # K
q_low, q_high = 5e3, 60e3  # search bounds (W/m)

# --- binary search for q' giving min margin ~ target ---
for _ in range(24):
    q_mid = 0.5*(q_low + q_high)
    res = run_channel(L,D_h,T_in,P,m_dot,q_mid,eps)
    min_margin = res["margin"].min()
    if min_margin >= target_margin:
        q_low = q_mid
    else:
        q_high = q_mid

q_star = q_low
res = run_channel(L,D_h,T_in,P,m_dot,q_star,eps)

print(f"q' ≈ {q_star/1e3:.1f} kW/m achieves min margin ≈ {res['margin'].min():.2f} K")
print(f"T_out = {res['T_c'][-1]:.2f} K, ΔP = {res['dP'][-1]/1e6:.3f} MPa")

# quick plot for README
base = Path(__file__).resolve().parents[2]
fig_dir = base / "figures"; fig_dir.mkdir(parents=True, exist_ok=True)
x = res["x"]
plt.figure()
plt.plot(x, res["margin"], label="ΔT_sat")
plt.axhline(0, ls="--")
plt.xlabel("x (m)"); plt.ylabel("ΔT_sat (K)")
plt.title(f"Margin with tuned q' = {q_star/1e3:.1f} kW/m")
plt.grid(True, linestyle=":")
plt.legend()
plt.savefig(fig_dir / "margin_tuned_qprime.png", dpi=200, bbox_inches="tight")
