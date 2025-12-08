import unittest
import numpy as np
from run_rodModels import Rod
from models.LinearRod import linear_model


class test_LinearRod(unittest.TestCase):
    def test_LinearRod1(self):
        """check whether the model works for a known case"""
        acc = 5                                                                     # order of accuracy of checks
        test_rod = Rod(E=1.0, L=1.0, I=1.0)                                         # properties of the test rod
        F1, F3, P, Q, N = 1.0, 1.0, 1.0, 1.0, 1.0                                   # loading condition of test case
        y_L, V_0, V_L, M_0, M_L = 0.29167, 0, 1, 0.5, 1                             # expected boundary values of test case


        test_rod.boundary_conditions(F1_mag=F1,F3_mag=F3,P_mag=P,Q_mag=Q,N_mag=N)
        linear_model(test_rod)

        self.assertAlmostEqual(test_rod.y[-1],y_L,places=acc)
        self.assertAlmostEqual(test_rod.f1[0],V_0,places=acc)
        self.assertAlmostEqual(test_rod.f1[-1],V_L,places=acc)
        self.assertAlmostEqual(test_rod.q2[0],M_0,places=acc)
        self.assertAlmostEqual(test_rod.q2[-1],M_L,places=acc)

    def test_LinearRod2(self):
        """check whether the model works for another known case"""
        acc=5                                                                       # order of accuracy of checks
        test_rod = Rod(E=15.24, L=1.54, I=0.76)                                     # properties of the test rod
        F1, F3, P, Q, N = 1.51, 0.0, 5.2, 23.5, 3.11                                # loading condition of test case
        y_L, V_0, V_L, M_0, M_L = 1.951005, 2.8746, 5.2, 17.282557, 23.5            # expected boundary values of test case

        test_rod.boundary_conditions(F1_mag=F1, F3_mag=F3, P_mag=P, Q_mag=Q, N_mag=N)
        linear_model(test_rod)

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
        linear_model(test_rod)

        self.assertAlmostEqual(test_rod.y[0],0)                             # fixed end should not move
        self.assertAlmostEqual(test_rod.f1[-1], V_L)
        self.assertAlmostEqual(test_rod.q2[-1], M_L)


if __name__ == '__main__':
    unittest.main()

