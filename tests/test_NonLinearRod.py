import unittest
from run_rodModels import Rod
from models.NonLinearRod import nonlinear_model
import numpy as np

class test_NonLinearRod(unittest.TestCase):
    def test_NonLinearRod1(self):
        """check whether the model works for a known case"""
        acc = 3                                                                     # order of accuracy of checks; accuracy reduced as solve_ivp is numerically unstable
        test_rod = Rod(E=1.0, L=1.0, I=1.0)                                         # properties of the test rod
        F1, F3, P, Q, N = 1.0, 0.0, 1.0, 1.0, 1.0                                   # loading condition of test case
        y_L, V_0, V_L, M_0, M_L = 0.389933, -0.9518812, 1, 0.97630532, 1            # expected boundary values of test case

        test_rod.boundary_conditions(F1_mag=F1,F3_mag=F3,P_mag=P,Q_mag=Q,N_mag=N)
        nonlinear_model(test_rod)

        self.assertAlmostEqual(test_rod.y[-1],y_L,places=acc)
        self.assertAlmostEqual(test_rod.f1[0],V_0,places=acc)
        self.assertAlmostEqual(test_rod.f1[-1],V_L,places=acc)
        self.assertAlmostEqual(test_rod.q2[0],M_0,places=acc)
        self.assertAlmostEqual(test_rod.q2[-1],M_L,places=acc)

    def test_NonLinearRod2(self):
        """check whether the model works for another known case"""
        acc=3                                                                       # order of accuracy of checks; accuracy reduced as solve_ivp is numerically unstable
        test_rod = Rod(E=15.24, L=1.54, I=0.76)                                     # properties of the test rod
        F1, F3, P, Q, N = 1.51, 0.0, 5.2, 23.5, 3.11                                # loading condition of test case
        y_L, V_0, V_L, M_0, M_L = 0.982325, -5.097810, 5.2, 27.036914, 23.5         # expected boundary values of test case

        test_rod.boundary_conditions(F1_mag=F1, F3_mag=F3, P_mag=P, Q_mag=Q, N_mag=N)
        nonlinear_model(test_rod)

        self.assertAlmostEqual(test_rod.y[-1], y_L, places=acc)
        self.assertAlmostEqual(test_rod.f1[0], V_0, places=acc)
        self.assertAlmostEqual(test_rod.f1[-1], V_L, places=acc)
        self.assertAlmostEqual(test_rod.q2[0], M_0, places=acc)
        self.assertAlmostEqual(test_rod.q2[-1], M_L, places=acc)

    def test_bc(self):
        """check whether the boundary conditions are satisfied after the solution for any random value"""
        rng = np.random.default_rng()
        F1, F3, P, Q, N = rng.random(5)
        V_L, M_L = P, Q
        E, L, I = rng.random(3)
        test_rod = Rod(E=E, L=L, I=I)
        test_rod.boundary_conditions(F1_mag=F1, F3_mag=F3, P_mag=P, Q_mag=Q, N_mag=N)
        nonlinear_model(test_rod)

        self.assertAlmostEqual(test_rod.y[0],0)                                 # fixed end should not move
        self.assertAlmostEqual(test_rod.f1[-1], V_L)
        self.assertAlmostEqual(test_rod.q2[-1], M_L)


if __name__ == '__main__':
    unittest.main()

