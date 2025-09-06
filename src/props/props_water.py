"""Water-like properties at ~10 MPa for educational use (constant values)."""
import numpy as np
def water_props_const(T=560.0, P=10e6):
    rho = 750.0        # kg/m^3
    mu = 1.5e-4        # Pa·s
    k = 0.5            # W/m-K
    cp = 5200.0        # J/kg-K
    Pr = cp * mu / k
    return dict(rho=rho, mu=mu, k=k, cp=cp, Pr=Pr)


# --- NEW: crude T-dependent fits near ~10 MPa (educational only) ---
def water_props_T(T, P=10e6):
    """
    Very rough, smooth functions for 10 MPa water-like properties.
    Valid-ish over ~520–800 K for educational modeling only.
    """
    # density: weakly decreasing with T
    rho = 780.0 - 0.07*(T - 560.0)          # kg/m^3
    rho = max(rho, 600.0)

    # viscosity: drop with T
    mu0 = 2.2e-4
    mu  = mu0 * np.exp(-0.012*(T - 560.0))  # Pa·s
    mu  = max(mu, 6.0e-5)

    # conductivity: slight fall with T
    k   = 0.58 - 1.2e-4*(T - 560.0)         # W/m-K
    k   = max(k, 0.35)

    # cp: modest increase with T
    cp  = 4800.0 + 1.8*(T - 560.0)          # J/kg-K

    Pr  = cp * mu / k
    return dict(rho=rho, mu=mu, k=k, cp=cp, Pr=Pr)