"""Reduced benchmark constants and exact numerical reference used in the revision."""

import casadi as ca

import numpy as np

DH = 162000.0

rho = 900.0

cp = 2000.0

U = 20.0

sigma = 0.076

k0 = 97200.0

Ea = 82500.0

R = 8.314

Na0 = 2.5

Nb0 = 2.5

Va0 = 0.00032

Vb0 = 0.00023

Tmax = 397.0

T = 350.0

Tc_min = 335.0

u_min = 0.0

u_max = 2.777777777777778e-08

Vmax = 0.00055

tf = 108000.0

k_T = 4.729015269131589e-08

cb_feed = 10869.565217391304

Kp = 5e-07

Ki = 3e-11

e_sp = 1.0

r_plus = 2.777777777777778e-08

r_minus = 8.333333333333334e-09

eta_s = 0.15

eta_e = 0.05

def Na_v(xa):
    return Na0 * (1 - xa)


def Nb_v(xa, V):
    return Nb0 * (V - Va0) / Vb0 - Na0 * xa


def rate(xa, V):
    nb = max(Nb_v(xa, V), 0.0)
    return k_T * Na_v(xa) * nb / V**2


def Tcf(xa, V):
    nb = max(Nb_v(xa, V), 0.0)
    return T + DH * min(Na_v(xa), nb) / (rho * cp * V)


def solve_ocp(N=100):
    dt = tf / N
    opti = ca.Opti()
    X = opti.variable(2, N + 1)
    Uc = opti.variable(1, N)

    def _Na(xa):
        return Na0 * (1 - xa)

    def _Nb(xa, V):
        return Nb0 * (V - Va0) / Vb0 - Na0 * xa

    def _rate(xa, V):
        return k_T * _Na(xa) * _Nb(xa, V) / V**2

    def _Tcf(xa, V):
        return T + DH * ca.fmin(_Na(xa), _Nb(xa, V)) / (rho * cp * V)

    def _Qrx(xa, V):
        return DH * _rate(xa, V) * V

    def _Qmax(V):
        return U * (2 * V / sigma + ca.pi * sigma**2) * (T - Tc_min)

    def _ode(x, u):
        return ca.vertcat(_rate(x[0], x[1]) * x[1] / Na0, u)

    opti.minimize(-X[0, N])
    opti.subject_to(X[0, 0] == 0)
    opti.subject_to(X[1, 0] == Va0)

    for k in range(N):
        xk = X[:, k]
        uk = Uc[:, k]
        k1 = _ode(xk, uk)
        k2 = _ode(xk + dt / 2 * k1, uk)
        k3 = _ode(xk + dt / 2 * k2, uk)
        k4 = _ode(xk + dt * k3, uk)
        opti.subject_to(X[:, k + 1] == xk + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4))

    opti.subject_to(opti.bounded(u_min, Uc, u_max))
    for k in range(N + 1):
        opti.subject_to(X[1, k] <= Vmax)
        opti.subject_to(X[1, k] >= Va0 * 0.99)
        opti.subject_to(_Tcf(X[0, k], X[1, k]) <= Tmax - 1e-3)
        opti.subject_to(_Qrx(X[0, k], X[1, k]) <= _Qmax(X[1, k]) - 1e-3)

    opti.set_initial(Uc, u_max * 0.3)
    for k in range(N + 1):
        opti.set_initial(X[0, k], k / N * 0.5)
        opti.set_initial(X[1, k], Va0 + (Vmax - Va0) * k / N)

    opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "max_iter": 5000, "tol": 1e-9})
    sol = opti.solve()
    return (
        np.linspace(0, tf, N + 1),
        sol.value(X[0, :]),
        sol.value(X[1, :]),
        sol.value(Uc).flatten(),
    )
