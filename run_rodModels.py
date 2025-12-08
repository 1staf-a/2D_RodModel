import matplotlib.pyplot as plt
import numpy as np

from models.LinearRod import linear_model
from models.NonLinearRod import nonlinear_model

class Rod: #2D cylindrical rod
    """
    A discrete 2D rod constrained to one end with a fixed support with distributed across its length and concentrated loads on its free end
    Used to simulate the deflection of a beam with linear and nonlinear theories of elasticity
    """
    def __init__(self,L=1.0,Nu=10,E=1.0,I=1.0,plt_arrows=False,modelName=""):
        """
        Initiate the rod object
        :param L: Length of the rod
        :param Nu: Number of nodes in the rod
        :param E: Young's modulus of the rod
        :param I: Moment of Inertia of the rod
        :param plt_arrows: flag for plotting shear loading arrows
        :param modelName: Name of elasticity model being used on the rod (Linear/Nonlinear)
        """
        self.model = modelName
        self.L = float(L)                                       # default rod length to 1[m]
        self.N = int(Nu)                                        # number of sections
        self.ds = float(L/(self.N-1))                           # minimum distance between cross-sections [m]

        self.s = np.linspace(0, L, self.N)                 # parametric length [m]

        self.I = float(I)                                       # moment of inertia of cross-sections [m^4]
        self.E = float(E)                                       # elastic modulus [N/m2]
        self.y = np.zeros([self.N])                             # deflection of cross-sections in [m]
        self.x = np.linspace(0, L, self.N)                 # x position of cross-sections in [m]

        self.k2 = np.zeros([self.N])                            # curvature of each cross-section

        self.r = np.hstack((self.x,self.y))                     # position vector of cross-section centreline [m,m]
        self.phi = np.zeros([self.N])                           # angles of the cross-sections centerline [rad]
        self.t = np.ones([self.N,2])                            # tangent unit vector for each cross-section
        self.n = np.ones([self.N,2])                            # normal unit vector of each cross-section

        self.F1 = 0                                             # constant distributed load [N/m] (Points +j)
        self.F3 = 0                                             # constant distributed load [N/m] (Points +i)
        self.P=0                                                # force at free end [N] (Points -j)
        self.Nf=0                                               # normal force at free end [N] (Points +i)
        self.Q=0                                                # moment at free end [Nm] (Points +k)

        self.f1 = np.zeros([self.N,1])                          # internal shear loading [N]
        self.f3 = np.zeros([self.N,1])                          # internal normal loading [N]
        self.q2 = np.zeros([self.N,1])                          # internal moment [Nm]

        self.plot_arrows=plt_arrows                             # flag for plotting quiver arrows

        self.eps=1e-6                                           # value for avoiding zero in denominators

    def plot_self(self):
        """
            Plot the rod object, shear loads and internal moments in a single subplot.
            :param self: inputs the rod object which has all the defined loads applied to it.
            :return:
            """
        fig = plt.figure()
        fig.set_size_inches(10, 10)
        ax1 = plt.subplot(3, 1, 1)

        ax1.grid()

        plt.title(self.model)
        max_def = np.amax(np.abs(self.y))
        plt.plot(self.x.T, self.y.T)
        if self.plot_arrows:
            plt.scatter(self.x, self.y)
            if self.model=="Linear":
                plt.quiver(self.x, self.y, 0*self.n[:, 0], self.f1/(np.max(np.abs(self.f1))+self.eps), scale=20, width=0.0045, color='k')
            else:
                plt.quiver(self.x, self.y, self.n[:, 0]*self.f1/(np.max(np.abs(self.f1))+self.eps), self.n[:, 1]*self.f1/(np.max(np.abs(self.f1))+self.eps), scale=20, width=0.0045, color='k')

        plt.tick_params('x', labelbottom=True)
        plt.ylabel('Y [m]')
        plt.ylim(max_def * -2, max_def * 2)
        plt.xlabel('X [m]')

        ax2 = plt.subplot(312)
        ax2.grid()
        plt.plot(self.s.T, self.f1.T)
        plt.tick_params('x', labelbottom=True)
        plt.ylabel('Internal Shear [N]')
        plt.xlabel('Rod Length [m]')

        ax3 = plt.subplot(313)
        ax3.grid()
        plt.plot(self.s.T, self.q2.T)
        plt.ylabel('Internal Moment [Nm]')
        plt.xlabel('Rod Length [m]')
        plt.show()

    def boundary_conditions(self,F1_mag=1.0,F3_mag=0.0,P_mag=0.0,Q_mag=0.0,N_mag=0.0):
        """
            Apply the boundary conditions to the rod
            :param self: provide a rod class object
            :param F1_mag: scalar multiple(N/m) for the distributed load being applied. (+j)-direction
            :param F3_mag: scalar multiple(N/m) for the distributed load being applied. (+i)-direction
            :param P_mag: concentrated load(N) being applied at the free end (-j)-direction
            :param Q_mag: concentrated moment(Nm) being applied at the free end (+k)-direction
            :param N_mag: concentrated load(N) being applied at the free end (-i)-direction
            :return:
            """

        self.F1 = F1_mag
        self.F3 = F3_mag
        self.P = P_mag
        self.Q = Q_mag
        self.Nf = N_mag

def compare_models(*rods):
    """
    Compare multiple rods objects using plots. Plots their deflections, shear loading and internal moments
    :param rods: rods objects
    :return:
    """
    legend = []
    plt.figure(0)
    plt.title('Deflections')
    plt.xlabel("X [m]")
    plt.ylabel("Y [m]")
    plt.figure(1)
    plt.title('Internal Shears')
    plt.xlabel("Rod Length [m]")
    plt.ylabel("Internal Shear [N]")
    plt.figure(2)
    plt.title('Internal Moments')
    plt.xlabel("Rod Length [m]")
    plt.ylabel("Internal Moment [Nm]")

    colors=["#4CAF50","#2196F3","#FFC107"]
    line_styles=["solid","dashed","dashdot"]
    i=0
    for c_rod in rods:
        legend.append(c_rod.model)
        plt.figure(0)
        plt.plot(c_rod.x, c_rod.y,label=c_rod.model,color=colors[i],linestyle=line_styles[i])
        if c_rod.plot_arrows:
            if c_rod.model == "Linear":
                plt.quiver(c_rod.x, c_rod.y, 0, 1*(c_rod.f1/(np.max(np.abs(c_rod.f1))+c_rod.eps)),
                           scale=21, width=0.0035,alpha=0.75,label=None,color=colors[i])
            else:
                plt.quiver(c_rod.x, c_rod.y, c_rod.n[:, 0]*(c_rod.f1/(np.max(np.abs(c_rod.f1)) + c_rod.eps)), c_rod.n[:, 1]*
                       (c_rod.f1/(np.max(np.abs(c_rod.f1))+c_rod.eps)), scale=21, width=0.0035,alpha=0.75,label=None,color=colors[i])
        plt.legend()
        plt.figure(1)
        plt.plot(c_rod.s, c_rod.f1,label=c_rod.model,color=colors[i],linestyle=line_styles[i])
        plt.legend()
        plt.figure(2)
        plt.plot(c_rod.s, c_rod.q2,label=c_rod.model,color=colors[i],linestyle=line_styles[i])
        plt.legend()
        i+=1

    plt.show()

def get_input_float(var):
    """
    prompts the user for a valid input in the form of a floating point
    :param var: name of the desired variable
    :return: user input
    """
    y = True
    fl_x= float()
    while y:
        fl_x = input(f"{var}:")
        try:
            fl_x = float(fl_x)
            y = False
        except:
            print("Wrong input, please try again.")

    return fl_x

def user_inputs(test_case=False):
    """
    function to get the loading conditions and properties of the rod from the user
    :param test_case: for running the test case for debugging (True if running test case)
    :return: rod properties, loadings and whether to plot the normal vectors and rods
    """
    if test_case:
        return [[1.0, 21, 1.0, 1.0], [0.0,0.0,10.0,0.0,0.0],True, True]

    moment = get_input_float("Moment[Nm] at free end")
    force1 = get_input_float("Force[N] 'P' at free end")
    force3 = get_input_float("Force[N] 'Nf' at free end")
    d_force1 = get_input_float("Distributed force[N/m] along 'F1'")
    d_force3 = get_input_float("Distributed force[N/m] along 'Nf'")

    Elasticity = get_input_float("Elasticity[N/m^2]")
    if Elasticity == 0:
        print("Elasticity must be non zero positive, defaulting to 1.")
        Elasticity = 1
    elif Elasticity < 0:
        print("Elasticity must be non-zero positive, defaulting to absolute.")
        Elasticity = abs(Elasticity)

    Inertia = get_input_float("Moment of inertia[m^4]")
    if Inertia == 0:
        print("Inertia must be non zero positive, defaulting to 1.")
        Inertia = 1
    elif Inertia < 0:
        print("Inertia must be non-zero positive, defaulting to absolute.")
        Inertia = abs(Inertia)

    Length = get_input_float("Length of rod[m]")
    if Length == 0:
        print("Length must be non zero positive, defaulting to 1.")
        Length = 1
    elif Length < 0:
        print("Length must be non-zero positive, defaulting to absolute.")
        Length = abs(Length)

    p_arrows = input("Plot arrows (y/n):")
    if p_arrows == 'y':plot_arrows = True
    else:plot_arrows = False

    p_rods = input("Plot Individual Rods (y/n):")
    if p_rods == 'y':
        plot_rods = True
    else:
        plot_rods = False

    # analytically calculate the max deflection at free end for the linear rod model
    y_max = (d_force1 * Length ** 4 / 8 + force1 * Length ** 3 / 3 + moment * Length ** 2 / 2) / (Elasticity * Inertia)

    # determine the number of segments from max deflection to an upper limit of 101
    max_segments = 1000
    min_segments = 21
    exp_coeff = 0.157
    Segments = np.min([min_segments*np.exp(exp_coeff*y_max),max_segments]).astype(int)

    properties = np.array([Length, Segments, Elasticity, Inertia], dtype=np.float64)
    loadings = np.array([d_force1, d_force3, force1, moment, force3], dtype=np.float64)

    return properties, loadings, plot_arrows, plot_rods

if __name__ == '__main__':
    Test_case = input("Run Test Case(y/n):")
    if Test_case == 'y': Run_test = True
    else: Run_test = False
    Properties, Loadings, Plot_arrows, Plot_rods = user_inputs(test_case=Run_test)

    rod = Rod(*Properties,Plot_arrows,modelName="Linear")
    rod2= Rod(*Properties,Plot_arrows,modelName="NonLinear")

    rod.boundary_conditions(*Loadings)
    linear_model(rod)

    rod2.boundary_conditions(*Loadings)
    nonlinear_model(rod2)

    if Plot_rods:
        rod.plot_self()
        rod2.plot_self()

    compare_models(rod,rod2)



