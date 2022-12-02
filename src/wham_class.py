#!/usr/bin/env python
import pickle
import numpy as np
from numpy import inf, nan

class Wham:
    """This is a class for doing WHAM and WHAM-D

    This is a class for doing the Weighted Histogram Analysis Method (WHAM) and its
    derivative through fluctuation theory (WHAM-D)

    Attributes:
        kbT (float): Thermal energy constant in current units
        rlow (float): Low distance cutoff in current units
        rhi (float): High distance cutoff in current units
        k (float): Force constant in current units
        nwindows (int): Number of WHAM Windows
        nbins (int): Number of WHAM histogram bins
        xc (array_like): Location of bias centers shape(nwindows,)
        xvals (array_like): Location of bin centers shape(nbins,)
        U (array_like): Potential energy of bias shape(nwindows,nbins)
        ni (array_like): Histogram counts shape(nwindows, nbins)
        cnt (array_like): Total histogram counts shape(nwindows,)
        hval (array_like): Normalized histogram shape(nwindows,nbins)
        center (array_like): Location of centers shape(nbins,)
        P (array_like): Potential energy estimate shape(nbins,)
        F (array_like): Free energy estimate shape(nwindows,)
        F_old (array_like): Old free energy estimate shape(nwindows,)
        isconverged (bool): Flag for whether converged
        eweight (bool): Flag for whether to do WHAM-D 
        whamcomplete (bool): Flag for whether WHAM has completed.



    """

    def __init__(self, xc, k, rlow, rhi, nwindows, nbins, tolerance = 1e-4,
                 kb=0.0019872041, T=298.15, eweight=False, bias=1):
        """Initializes the class and builds attributes

        This class takes built in data and builds the class attributes.

        Args:
            xc (array_like): Positions of the bin centers
            k (float): Force constant in correct units
            rlow (float): Lowest r position
            rhi (float): Highest r position
            nwindows (int): Number of WHAM windows
            nbins (int): Number of histogram bins
            tolerance (float): Tolerance for the WHAM calculation [default=1e-4]
            kb (float): Boltzmann's constant in current, consistent, units [default=0.0019872041]
            T (float): Temperature in current, consistent, units [default=298.15]
            eweight (bool): Flag for whether to use fluctuation theory [default=False]
        
        """
        self.kbT = kb*T
        self.tolerance = 1E-4
        self.nwindows = nwindows
        self.nbins = nbins
        self.rlow, self.rhi = rlow, rhi
        self.xc = xc
        self.k = k
        self.bias = bias

        # Derived attributes
        self.xvals = np.linspace(self.rlow, self.rhi, self.nbins)


        # Main Arrays
        self.U = []
        self.ni = []
        self.cnt = []
        self.hval = []
        self.center = []
        self.P = []

        # Derived values
        self.F_old = np.zeros(self.nwindows)
        self.F = self.F_old

        # Flags
        self.isconverged = False
        self.eweight = eweight
        self.whamcomplete = False

        # Do initial setup
        self.Build_Data()

        return

    def Build_Data(self):
        def choose_bias(bias):
            if bias == 1:
                def calc_bias(x,k):
                    # This function calculates the bias potential value 
                    return k*x**2.
                return calc_bias
            elif bias == 2:
                def calc_bias(x,k):
                    # This function calculates the bias potential value
                    return 0.5*k*x**2.
                return calc_bias
            else:
                exit("improper bias")

        calc_bias = choose_bias(self.bias)

        for window in range(self.nwindows):
            if window%10 == 0: print(window)
            data, en, den, dhist = [], {}, {}, {}

            try:
                data = pickle.load(open("wham_pckl/window_%d.pckl"%window,'rb'))
            except:
                raise OSError("Error: Need to generate pckl files first!")

            if self.eweight == True:
                try:
                    en = pickle.load(open("wham_pckl/ener_%d.pckl"%window,'rb'))
                except:
                    raise OSError("Error: Need to generate energy pckl first")
            
            hist, bins = np.histogram(data, bins=self.nbins, range=(self.rlow,self.rhi), density=False)

            if self.eweight == True:
                for key in en:
                    en[key] = en[key]
                    den[key] = den[key]-np.average(en[key])
                    dhist[key], dbins = np.histogram(data, bins=self.nbins, range=(self.rlow,self.rhi),
                                                     density=False, weights=-den[key])
                    dhval[key].append(dhist[key]/np.sum(hist))
            
            # Store in the class
            self.center = (bins[:-1] + bins[1:]) / 2
            self.ni.append(hist)
            self.hval.append(hist/np.sum(hist))
            self.cnt.append(np.sum(hist))
            self.dx = np.subtract(self.xc[window], self.center)
            self.U.append(calc_bias(self.dx,self.k))

        self.U = np.array(self.U)

    def Do_WHAM(self,maxiter=10000,tol=1e-6):
        """
        """
        iteration = 0
        while self.isconverged == False:
            self.Wham_Iteration()
            iteration += 1
            if iteration > maxiter:
                exit("Error: Reached too many iterations. Increase maxiter.")
        return

    def Wham_Iteration(self):
        """Does one WHAM iteration
        
        Does a regular WHAM iteration

        Args:
            F (array_like): Free energy estimate

        """
        self.Update_P()
        self.Update_F()

        # Calculate the error
        Ferr = np.sum(np.abs(np.subtract(self.F,self.F_old)))
        self.F_old = self.F
        print("Error: %s" % Ferr)
        
        #Check if converged
        if Ferr < self.tolerance:
            self.isconverged = True

        return

    def Update_P(self):
        """Function to calculate the probability distribution

        This function calculates the probability distribution, given the free energy
        estimate, F, and various other quantities.

        """
        Ut = self.U.T
        Fm = self.F - np.min(self.F)

        numerator = np.sum(self.ni,axis=0)

        # Build Denominator
        inside_exponent = np.divide(-np.subtract(Ut, Fm), self.kbT)
        pre_sum = np.multiply(self.cnt, np.exp(inside_exponent))
        denominator = np.sum(pre_sum,axis=1)

        self.P = numerator/denominator

    def Update_F(self):
        """Calculates the WHAM free energy estimate

        This function calculates the WHAM free energy estimate, F as 
        part of the self consistent solution. 

        """
        inside_sum = np.multiply(self.P, np.exp(-self.U/self.kbT))
        self.F = -self.kbT*np.log(np.sum(inside_sum,axis=1))
        
        # Remove infinity
        self.F[self.F==inf] = 0.0
    
    def Calc_PMF(self):
        pmf = -self.kbT*np.log(P)-np.min(-self.kbT*np.log(self.P))
        return pmf
    
    def Plot_PMF(self):
        import matplotlib.pyplot as plt
        plt.figure(dpi=300,figsize=(3,3))
        plt.plot(self.xvals,self.Report_PMF())
        plt.show()
    



