import numpy as np
from scipy.integrate import cumulative_trapezoid

def linear_model(rod_i):
    """
    Solves the linear elastic rod model for a perpendicular distributed(force) or concentrated load(at free end: force or moment)
    :param rod_i: inputs the rod object which has all the defined loads applied to it.
    :return: the deformed rod object.
    """
    rod_i.f1 = -(rod_i.F1 * (rod_i.L - rod_i.s) - rod_i.P)
    rod_i.f3 = -(rod_i.F3 * (rod_i.L - rod_i.s) - rod_i.Nf)
    c2 = (rod_i.Q - rod_i.P * rod_i.L + 0.5 * rod_i.F1 * rod_i.L ** 2)
    rod_i.q2 = rod_i.P * rod_i.s - rod_i.F1 * rod_i.L * rod_i.s + 0.5 * rod_i.F1 * rod_i.s ** 2 + c2
    rod_i.k2 = rod_i.q2 / (rod_i.E * rod_i.I)
    rod_i.x = rod_i.s
    rod_i.y = (rod_i.F1*rod_i.s**4/24+(rod_i.P-rod_i.F1*rod_i.L)*rod_i.s**3/6+(rod_i.Q-rod_i.P*rod_i.L+0.5*rod_i.F1*rod_i.L**2)*rod_i.s**2*0.5)/(rod_i.I*rod_i.E)

    rod_i.r = np.array((rod_i.x, rod_i.y)).T
    rod_i.phi[1:] = np.atan((rod_i.y[1:]-rod_i.y[:-1])/(rod_i.x[1:]-rod_i.x[:-1]))

    rod_i.t = np.array([np.cos(rod_i.phi), np.sin(rod_i.phi)]).T
    rod_i.n = np.array([-np.sin(rod_i.phi), np.cos(rod_i.phi)]).T

    # alternative method to compute deflections and rod positions; constrains the total length of the rod!
    # rod_i.phi[1:] = cumulative_trapezoid(rod_i.k2, rod_i.s)
    # rod_i.x[1:] = cumulative_trapezoid(np.cos(rod_i.phi), rod_i.s)
    # rod_i.y[1:] = cumulative_trapezoid(np.sin(rod_i.phi), rod_i.s)
    pass
