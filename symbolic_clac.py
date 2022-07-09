import sympy as sp

M, gamma, P_P0, rho_rho0, T_T0 = sp.symbols("M, gamma, P_P0, rho_rho0, T_T0")

eqn_T_T0 = ((1 + ((gamma - 1)/2) * (M ** 2)) ** (-1)) - T_T0
sol_T_T0 = sp.solve(eqn_T_T0, M)[0]
print(f"M(T_T0) = {sol_T_T0}")

eqn_P_P0 = ((1 + ((gamma - 1)/2) * (M ** 2)) ** (-gamma/(gamma-1))) - P_P0
sol_P_P0 = sp.solve(eqn_P_P0, M)[0]
print(f"M(P_P0) = {sol_P_P0}")

eqn_rho_rho0 = ((1 + ((gamma - 1)/2) * (M ** 2)) ** (-1/(gamma-1))) - rho_rho0
sol_rho_rho0 = sp.solve(eqn_rho_rho0, M)[0]
print(f"M(rho_rho0) = {sol_rho_rho0}")
