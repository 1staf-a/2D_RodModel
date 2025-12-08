# 2D_RodModel
A rod model simulating a 2D cantilever beam under user-defined loading conditions and compares the outcomes for a linear and nonlinear beam deflection model.

Python model of the publication: "AN INTERACTIVE GUI FOR EXPLORING CANTILEVER BENDING THROUGH NONLINEAR ROD THEORY", by Muhammad Hassaan Ahmed, Utkarsh Srivastava, Ranjan Das and Sachin Goyal, 2025, ASME IDETC.CIE2025.
Adapted from [repository](https://github.com/utkarshdrh/Static-Rod-Model-Codes-and-SImulink-Model) for Python.

# Instructions for running code
Execute /run_rodModels.py

User will be prompted for running test case. Test case runs sample simulation with predefined rod properties and simulation settings.

If test case is not run, user will be prompted to enter loading conditions and rod properties. Loads are defined as shown in figure:

<p align="center">
  <img width="826" height="257" alt="Rod" src="https://github.com/user-attachments/assets/156b4b56-95ae-48e5-b911-a4866d01c65a" />
</p>

Where P & Nf are concentrated applied forces on the free end, M is the concentrated applied moment at the free end. F1 & F3 are the distributed loads.
E is the Young's modulus of the rod, I is the moment of inertia of its cross section and L is the initial length of the rod.
Inputting invalid values(<=0) for E, I & L will resort to using the default in-built values.

User will be prompted to draw the shear arrows which will plot the shear loads on the deflection curves.
User will also be prompted to add another visualization which draws the linear and non-linear models separately with their respective internal loads.

Once the inputs have been defined the code will run for desired conditions and show the plots for the input conditions. 


## Large Deflection Conditions

For some input conditons which result in large deflections in the linear model, the curve of the non-linear model may not be visible. This not a bug, it is a consequence of the non-linear model conserving the total length of the rod while the linear model does not conserve the length of the rod.
