import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import cumulative_trapezoid

def nonlinear_model(rod_nl):
    """
    Solves the nonlinear rod model for a perpendicular distributed(force) or concentrated load(at free end: force or moment)
    :param rod_nl: inputs the rod object which has all the defined loads applied to it.
    :return: the deformed rod object.
    """
    def ode(s, K):
        # K = [k, k1, k2] = [k, k', k'']
        eps1 = 1e-8
        k = K[0]
        k1 = K[1]
        k2 = K[2]
        # ensure that it isn't dividing by zero!
        if k==0:
            k3 = 1 / (rod_nl.E * rod_nl.I * (k + eps1)) * (
                    rod_nl.E * rod_nl.I * k2 * k1 - rod_nl.E * rod_nl.I * k1 * k ** 3 - rod_nl.F1 * k1 - rod_nl.F3 * k ** 2)
        else:
            k3 = 1 / (rod_nl.E * rod_nl.I * k) * (
                    rod_nl.E * rod_nl.I * k2 * k1 - rod_nl.E * rod_nl.I * k1 * k ** 3 - rod_nl.F1 * k1 - rod_nl.F3 * k ** 2)

        return np.array([k1, k2, k3])

    K00= rod_nl.Q / rod_nl.E / rod_nl.I
    K01= -rod_nl.P / rod_nl.E / rod_nl.I
    K02= rod_nl.Q * rod_nl.Nf / (rod_nl.E * rod_nl.I) ** 2 + rod_nl.F1 / (rod_nl.E * rod_nl.I)
    K0=[K00,K01,K02]

    odesol = solve_ivp(ode, (0, rod_nl.L), K0, t_eval=rod_nl.s)
    rod_nl.k2 = np.flip(odesol.y[0, :])
    rod_nl.phi[1:]= cumulative_trapezoid(rod_nl.k2, rod_nl.s)
    rod_nl.x[1:] = cumulative_trapezoid(np.cos(rod_nl.phi), rod_nl.s)
    rod_nl.y[1:] = cumulative_trapezoid(np.sin(rod_nl.phi), rod_nl.s)

    rod_nl.q2 = rod_nl.E * rod_nl.I * rod_nl.k2
    rod_nl.f1 = -rod_nl.E * rod_nl.I * np.flip(odesol.y[1, :])

    eps = 1e-8
    rod_nl.k2[np.abs(rod_nl.k2)<=eps]=eps
    rod_nl.f3 = -(rod_nl.F1 + rod_nl.E * rod_nl.I * np.flip(odesol.y[2, :])) / (rod_nl.k2 + eps)

    rod_nl.r = np.array((rod_nl.x, rod_nl.y)).T
    rod_nl.t = np.array([np.cos(rod_nl.phi), np.sin(rod_nl.phi)]).T
    rod_nl.n = np.array([-np.sin(rod_nl.phi), np.cos(rod_nl.phi)]).T
    pass

def mat_nonlinear_model(rod_nl):
    """
    alternative nonlinear rod model for a perpendicular distributed(force) or concentrated load(at free end: force or moment)
    works by using the prescribed odes rather than simplifying to kappa form
    :param rod_nl: inputs the rod object which has all the defined loads applied to it.
    :return: the deformed rod object.
    """
    def odes(s, x):
        # E = [k2, f1, f3]
        e1 = -x[1] / rod_nl.E / rod_nl.I
        e2 = -x[0] * x[2] - rod_nl.F1
        e3 = x[0] * x[1] - rod_nl.F3
        return np.array([e1, e2, e3])
    E0=[rod_nl.Q/rod_nl.E/rod_nl.I,rod_nl.P,rod_nl.Nf]
    odesol = solve_ivp(odes, (0, rod_nl.L), E0, t_eval=rod_nl.s)
    rod_nl.k2 = np.flip(odesol.y[0, :])
    rod_nl.phi[1:] = cumulative_trapezoid(rod_nl.k2, rod_nl.s)
    rod_nl.x[1:] = cumulative_trapezoid(np.cos(rod_nl.phi), rod_nl.s)
    rod_nl.y[1:] = cumulative_trapezoid(np.sin(rod_nl.phi), rod_nl.s)

    rod_nl.q2 = rod_nl.E * rod_nl.I * rod_nl.k2
    rod_nl.f1 = np.flip(odesol.y[1, :])
    rod_nl.f3 = np.flip(odesol.y[2, :])

    rod_nl.r = np.array((rod_nl.x, rod_nl.y)).T
    rod_nl.t = np.array([np.cos(rod_nl.phi), np.sin(rod_nl.phi)]).T
    rod_nl.n = np.array([-np.sin(rod_nl.phi), np.cos(rod_nl.phi)]).T

    pass