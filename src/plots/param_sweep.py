import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Import the model parts from your film-T script
from src.props.props_water import water_props_T

# Reuse the helpers from axial_steady (copy here to keep this file standalone)
def h_from_props(rho, mu, k, cp, D_h, m_dot, A_flow):
    v  = m_dot / (rho * A_flow)
    Re = rho * v * D_h / mu
    Pr = cp * mu / k
    Nu = 0.023 * (Re**0.8) * (Pr**0.4)
    h  = Nu * k / D_h
    return h, Re, Nu, v

def f_blasius(Re):
    return 0.3164 / (Re**0.25)

def tsat_water(P_pa: float) -> float:
    try:
        import CoolProp.CoolProp as CP
        return CP.PropsSI("T", "P", P_pa, "Q", 0, "Water")
    except Exception:
        return 584.0  # placeholder near 10 MPa

def run_channel(L, D_h, T_in, P, m_dot, q_lin, N=200):
    x = np.linspace(0, L, N+1); dx = x[1]-x[0]
    P_wet = np.pi * D_h
    A_flow = np.pi * (D_h/2)**2
    q_pp = q_lin / P_wet

    T_c = np.zeros_like(x); T_c[0] = T_in
    T_w = np.zeros_like(x)
    h   = np.zeros_like(x)
    Re  = np.zeros_like(x)
    dP  = np.zeros_like(x)

    props0 = water_props_T(T_in, P)
    h_guess, Re0, _, v0 = h_from_props(
        rho=props0["rho"], mu=props0["mu"], k=props0["k"], cp=props0["cp"],
        D_h=D_h, m_dot=m_dot, A_flow=A_flow
    )
    T_w[0] = T_c[0] + q_pp / h_guess
    h[0], Re[0] = h_guess, Re0

    for i in range(N):
        props_i = water_props_T(0.5*(T_c[i]+T_w[i]), P)
        T_c[i+1] = T_c[i] + (q_lin/(m_dot*props_i["cp"])) * dx

        T_c_loc = T_c[i+1]
        T_w_loc = T_c_loc + q_pp / max(h[i], 1e-3)
        pr = props_i
        for _ in range(3):
            T_f = 0.5*(T_c_loc + T_w_loc)
            pr  = water_props_T(T_f, P)
            h_loc, Re_loc, _, v_loc = h_from_props(
                rho=pr["rho"], mu=pr["mu"], k=pr["k"], cp=pr["cp"],
                D_h=D_h, m_dot=m_dot, A_flow=A_flow
            )
            T_w_loc = T_c_loc + q_pp / h_loc

        T_w[i+1] = T_w_loc
        h[i+1]   = h_loc
        Re[i+1]  = Re_loc
        f = f_blasius(Re_loc)
        dPdx = f * (pr["rho"] * v_loc**2 / 2.0) / D_h
        dP[i+1] = dP[i] + dPdx * dx

    return dict(x=x, T_c=T_c, T_w=T_w, h=h, Re=Re, dP=dP)

# ---------------- Sweep setup ----------------
L=6.0; D_h=0.012; T_in=560.0; P=10e6
m_dots = [1.0, 1.5, 2.0, 2.5, 3.0]                # kg/s
q_lines = [60e3, 150e3, 250e3, 350e3]             # W/m

records = []
for m_dot in m_dots:
    for q_lin in q_lines:
        res = run_channel(L, D_h, T_in, P, m_dot, q_lin)
        Tmax_w = res["T_w"].max()
        Tout_c = res["T_c"][-1]
        dP_tot = res["dP"][-1]
        records.append({
            "m_dot_kg_s": m_dot,
            "q_lin_W_m": q_lin,
            "T_w_max_K": Tmax_w,
            "T_c_out_K": Tout_c,
            "DeltaP_total_MPa": dP_tot/1e6
        })

df = pd.DataFrame(records)

# ---------------- Save & plot ----------------
base = Path(__file__).resolve().parents[2]
data_dir = base / "data"; data_dir.mkdir(exist_ok=True, parents=True)
fig_dir  = base / "figures"; fig_dir.mkdir(exist_ok=True, parents=True)

df.to_csv(data_dir / "sweep_mdots_qline.csv", index=False)

# Plot 1: Max wall T vs m_dot for each q'
plt.figure()
for q_lin in q_lines:
    dfi = df[df["q_lin_W_m"]==q_lin].sort_values("m_dot_kg_s")
    plt.plot(dfi["m_dot_kg_s"], dfi["T_w_max_K"], marker="o", label=f"q'={q_lin/1e3:.0f} kW/m")
plt.xlabel("m_dot (kg/s)")
plt.ylabel("Max wall temperature T_w,max (K)")
plt.title("Effect of ṁ on Max Wall Temperature")
plt.grid(True, linestyle=":")
plt.legend()
plt.savefig(fig_dir / "sweep_Twmax_vs_mdot.png", dpi=200, bbox_inches="tight")

# Plot 2: Total ΔP vs m_dot for each q'
plt.figure()
for q_lin in q_lines:
    dfi = df[df["q_lin_W_m"]==q_lin].sort_values("m_dot_kg_s")
    plt.plot(dfi["m_dot_kg_s"], dfi["DeltaP_total_MPa"], marker="o", label=f"q'={q_lin/1e3:.0f} kW/m")
plt.xlabel("m_dot (kg/s)")
plt.ylabel("Total ΔP (MPa)")
plt.title("Effect of ṁ on Pressure Drop")
plt.grid(True, linestyle=":")
plt.legend()
plt.savefig(fig_dir / "sweep_dP_vs_mdot.png", dpi=200, bbox_inches="tight")

print("Saved:",
      data_dir / "sweep_mdots_qline.csv",
      fig_dir / "sweep_Twmax_vs_mdot.png",
      fig_dir / "sweep_dP_vs_mdot.png", sep="\n")
