#!/usr/bin/env python
import pickle, warnings
import numpy as np
from numpy import inf, nan
import matplotlib.pyplot as plt

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
        U (array_like): Potential energy of bias shape(nwindows,nbins)
        ni (array_like): Histogram counts shape(nwindows, nbins)
        cnt (array_like): Total histogram counts shape(nwindows,)
        hval (array_like): Normalized histogram shape(nwindows,nbins)
        center (array_like): Location of centers shape(nbins,)
        P (array_like): Potential energy estimate shape(nbins,)
        F (array_like): Free energy estimate shape(nwindows,)
        F_old (array_like): Old free energy estimate shape(nwindows,)
        isconverged (bool): Flag for whether converged
        d_converged (bool): Flag for whether derivative converged
        eweight (bool): Flag for whether to do WHAM-D 
        whamcomplete (bool): Flag for whether WHAM has completed.
        iter_output (bool): Flag for printing iteration output to the screen



    """

    def __init__(self, xc, k, rlow, rhi, nwindows, nbins, tolerance = 1e-4,
                 kb=0.0019872041, T=298.15, eweight=False, bias=1, iter_output=True):
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
            iter_output (bool): Falg for whether to print iteration output [default=True]
        
        """
        self.kbT = kb*T
        self.tolerance = 1E-4
        self.nwindows = nwindows
        self.nbins = nbins
        self.rlow, self.rhi = rlow, rhi
        self.xc = xc
        self.k = k
        self.bias = bias


        # Main Arrays
        self.U = []
        self.ni = []
        self.cnt = []
        self.hval = []
        self.center = []
        self.P = []
        self.dP = []

        # Derived values
        self.F_old = np.zeros(self.nwindows)
        self.F = self.F_old
        self.dF_old = np.zeros(self.nwindows)
        self.dF = self.dF_old

        # Flags
        self.isconverged = False
        self.d_converged = False
        self.eweight = eweight
        self.whamcomplete = False
        self.iter_output = iter_output

        # Do initial setup
        self.Build_Data()

        return

    def Build_Data(self):
        """ This reads in data and builds histograms

        This program reads in the collective variable data, and the
        energies (if the derivative WHAM-D calculation is turned on) and then
        it histograms the data, and generates the structures needed to run a WHAM (and WHAM-D)
        calculation.

        Raises:
            OSError: Something went wrong with the pickle files read in.
            ValueError: Wrong input for the bias potential was chosen. 

        """
        def choose_bias(bias):
            """Sets the bias used for the collective variable"""
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
                raise ValueError("Error: Improper bias option selected.")

        calc_bias = choose_bias(self.bias)

        self.dhval = {} # initialize derivative histogram
        
        # Loop over windows
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
            
            # Histogram the data over the selected range. Do not normalize
            hist, bins = np.histogram(data, bins=self.nbins, range=(self.rlow,self.rhi), density=False)

            # Do derivative weighting, if necessary
            if self.eweight == True:
                for key in en:
                    if key not in self.dhval: 
                        self.dhval[key] = [] # Initialize to key if not there. 
                    den[key] = en[key]-np.average(en[key])
                    dhist[key], dbins = np.histogram(data, bins=self.nbins, range=(self.rlow,self.rhi),
                                                     density=False, weights=-den[key])
                    self.dhval[key].append(dhist[key]/np.sum(hist))
            
            # Store in the class
            self.center = (bins[:-1] + bins[1:]) / 2
            self.ni.append(hist)
            self.hval.append(hist/np.sum(hist))
            self.cnt.append(np.sum(hist))
            self.dx = np.subtract(self.xc[window], self.center)
            self.U.append(calc_bias(self.dx,self.k))

        # Make U a numpy array, so that it is transposable (Same for dhval)
        self.U = np.array(self.U)
        if self.eweight == True:
            for key in self.dhval:
                self.dhval[key] = np.array(self.dhval[key])

    def Do_WHAM(self,maxiter=10000):
        """Function to do WHAM

        This function does WHAM by repeatedly calling the Wham_Iteration function, which
        works by updating P and F stored within the class. The function stops if the maxiter
        parameter is reached.

        Args:
            maxiter (int): Maximum number of WHAM iterations

        Raises:
            RuntimeWarning: WHAM didn't converge within maxiter steps to specified tolerance.

        """
        iteration = 0
        while self.isconverged == False:
            self.Wham_Iteration()
            iteration += 1
            if iteration > maxiter:
                warnings.warn("Warning: Reached too many iterations. Increase maxiter.",RuntimeWarning)
                break
    
    def Do_WHAM_D(self, key, maxiter=10000):
        """Function to do WHAM-D

        This function calculates the derivative of WHAM by repeatedly calling the Wham_D_Iteration 
        function, which works by updating dP and dF stored within the class. The function stops if
        the maxiter parameter is reached. This must be done after Do_WHAM has completed.

        Args:
            key (string): Energy key name
            maxiter (int): Maximum number of WHAM-D iterations

        Raises:
            RuntimeError: WHAM-D didn't converge within maxiter steps to specified tolerance.
            AssertionError: WHAM wasn't converged prior to calculating the derivative.
    
        """
        
        # Raises the AssertionError if WHAM wasn't converged first.
        self._Test_Convergant()

        iteration = 0
        self.d_converged = False
        while (self.d_converged == False):
            self.Wham_D_Iteration(key)
            iteration += 1

            if iteration > maxiter:
                warnings.warn("Warning: Too many iterations in derivative", RuntimeWarning)
                break

    def Wham_Iteration(self):
        """Does one WHAM iteration
        
        Does a regular WHAM iteration, by updating P, then F.

        """
        self.Update_P()
        self.Update_F()

        # Calculate the error
        Ferr = np.sum(np.abs(np.subtract(self.F, self.F_old)))
        self.F_old = self.F
        if self.iter_output: print("\r Error: %s" % Ferr)
        
        # Check if converged
        if Ferr < self.tolerance:
            self.isconverged = True

    
    def Wham_D_Iteration(self, key):
        """Does one WHAM-D iteration

        Does a regular WHAM-D iteration, by updating dP, then dF.

        """

        # Update the derivatives
        self.Update_dP(key)
        self.Update_dF()


        # Calculate the Error
        dFerr = np.sum(np.abs(np.subtract(self.dF,self.dF_old)))
        self.dF_old = self.dF

        print("Error: %s " % dFerr)
        print("dF: %s" % np.sum(self.dF))

        # Test Convergence
        self.d_converged = False
        if dFerr < self.tolerance:
            self.d_converged = True

    def Update_P(self):
        """Function to calculate the probability distribution

        This function calculates the probability distribution, given the free energy
        estimate, F, and various other quantities.

        """
        Ut = self.U.T
        Fm = self.F - np.min(self.F)

        numerator = np.sum(self.ni, axis=0)

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

    def Update_dP(self, key):
        """Calculates the WHAM-D derivative of probability distribution

        This function calculates the WHAM-D derivative of the probability distribution,
        dP, as part of a self consistent solution.

        Args:
            key (str): Energy type for the iteration

        """
        Ut = self.U.T
        Fm = self.F - np.min(self.F)
        dHT = self.dhval[key].T

        # Build Denominator
        inside_sum = np.multiply(self.cnt, np.exp(-np.divide(np.subtract(Ut,Fm),self.kbT)))
        denominator = np.sum(inside_sum,axis=1)

        # Build Term 1
        term1 = np.sum(np.multiply(self.cnt,-dHT),axis=1)

        # Build Term 2
        exponent = np.exp(-np.divide(np.subtract(Ut,Fm),self.kbT))
        weight_exp = np.multiply(self.cnt,exponent)
        prefact_exp = np.multiply(weight_exp,(self.dF/self.kbT+Fm-Ut))
        term2 = self.P*np.sum(prefact_exp,axis=1)

        self.dP = (term1 - term2)/denominator
        
    def Update_dF(self):
        """Calculates the WHAM-D derivative of free energy

        This function calculates the WHAM-D derivative of the free energy,
        dP, as part of a self consistent solution.

        """

        exponent = np.exp(-self.U/self.kbT)
        # Build Term 1
        term1 = -self.kbT*self.F
        
        # Build Term 2
        inside_sum = np.multiply(np.subtract(self.dP, np.multiply(self.U, self.P)), exponent)

        term2_num = self.kbT*np.sum(inside_sum, axis=1)
        term2_den = np.sum(np.multiply(exponent, self.P), axis=1)

        self.dF = (term1 - term2_num/term2_den)

        print(np.sum(term1),np.sum(term2_num/term2_den))
        self.dF[self.dF==inf] = 0.0

    def Do_Alt_WHAM_D(self, key, maxiter=10000):
        """Function to do WHAM-D

        This function calculates the derivative of WHAM by repeatedly calling the Wham_D_Iteration 
        function, which works by updating dP and dF stored within the class. The function stops if
        the maxiter parameter is reached. This must be done after Do_WHAM has completed.

        Args:
            key (string): Energy key name
            maxiter (int): Maximum number of WHAM-D iterations

        Raises:
            RuntimeError: WHAM-D didn't converge within maxiter steps to specified tolerance.
            AssertionError: WHAM wasn't converged prior to calculating the derivative.
    
        """
        
        # Raises the AssertionError if WHAM wasn't converged first.
        self._Test_Convergant()

        iteration = 0
        self.d_converged = False
        while (self.d_converged == False):
            self.Alt_Wham_D_Iteration(key)
            iteration += 1

            if iteration > maxiter:
                warnings.warn("Warning: Too many iterations in derivative", RuntimeWarning)
                break

    def Alt_Wham_D_Iteration(self, key):
        """Does one WHAM-D iteration

        Does a regular WHAM-D iteration, by updating dP, then dF.

        """

        # Update the derivatives
        self.Update_Alt_dP(key)
        self.Update_Alt_dF()


        # Calculate the Error
        dFerr = np.sum(np.abs(np.subtract(self.dF,self.dF_old)))
        self.dF_old = self.dF

        print("Error: %s " % dFerr)
        print("dF: %s" % np.sum(self.dF))

        # Test Convergence
        self.d_converged = False
        if dFerr < self.tolerance*1e-4:
            self.d_converged = True

    def Update_Alt_dP(self, key):
        """Calculates the WHAM-D derivative of probability distribution

        This function calculates the WHAM-D derivative of the probability distribution,
        dP, as part of a self consistent solution.

        Args:
            key (str): Energy type for the iteration

        """
        Ut = self.U.T
        Fm = self.F - np.min(self.F)
        dHT = self.dhval[key].T

        # Build Denominator
        inside_sum = np.multiply(self.cnt, np.exp(-np.divide(np.subtract(Ut,Fm),self.kbT)))
        denominator = np.sum(inside_sum,axis=1)

        # Build Term 1
        #term1 = np.sum(np.multiply(self.cnt,-dHT),axis=1)
        exponent = np.exp(-np.divide(np.subtract(Ut,Fm),self.kbT))
        pre_exp = np.subtract(Ut,Fm)
        pdiff = np.subtract(dHT,np.multiply(pre_exp,exponent))
        inside_t1_sum = np.multiply(self.cnt,pdiff)
        term1 = np.sum(inside_t1_sum,axis=1)/denominator

        # Build Term 2
        #exponent = np.exp(-np.divide(np.subtract(Ut,Fm),self.kbT))
        #weight_exp = np.multiply(self.cnt,exponent)
        #prefact_exp = np.multiply(weight_exp,(self.dF/self.kbT+Fm-Ut))
        #term2 = self.P*np.sum(prefact_exp,axis=1)
        term2=self.P*np.sum(self.cnt*exponent*self.dF,axis=1)/denominator
        print(np.sum(term1),np.sum(term2))
        self.term1 = term1
        self.term2 = term2

        self.dP = -(term1 + term2)
        
    def Update_Alt_dF(self):
        """Calculates the WHAM-D derivative of free energy

        This function calculates the WHAM-D derivative of the free energy,
        dP, as part of a self consistent solution.

        """

        exponent = np.exp(-self.U/self.kbT)
        # Build Term 1
        #term1 = -self.kbT*self.F
        
        # Build Term 2
        #inside_sum = np.multiply(np.subtract(self.dP, np.multiply(self.U, self.P)), exponent)

        #term2_num = self.kbT*np.sum(inside_sum, axis=1)
        #term2_den = np.sum(np.multiply(exponent, self.P), axis=1)

        #self.dF = (term1 - term2_num/term2_den)

        #print(np.sum(term1),np.sum(term2_num/term2_den))
        #self.dF[self.dF==inf] = 0.0

        term1 = -self.kbT*self.F
        inside_brace=np.subtract(np.multiply(self.P,self.U),self.dP)
        term2_num = self.kbT*np.sum(np.multiply(exponent,inside_brace),axis=1)
        term2_den = np.sum(exponent*self.P,axis=1)

        self.dF = term1 + term2_num/term2_den
        self.dF[self.dF==inf] = 0.0


    
    def _Test_Convergant(self):
        """Simple functon that tests convergence

        Raises:
            AssertionError: The calculation is not converged.

        """
        try:
            assert self.isconverged == True
        except:
            raise AssertionError("Error: WHAM-D must follow a successful WHAM calculation")



    def Calc_PMF(self):
        """Calculates the potential of mean force

        This calculates the potential of mean force from the self consistent
        solution to WHAM, which must have converged.

        Raises:
            AssertionError: Calculation wasn't converged before trying to plot

        """
        self._Test_Convergant()
        pmf = -self.kbT*np.log(self.P)-np.min(-self.kbT*np.log(self.P))
        return pmf
    
    def Plot_PMF(self):
        """Plots PMF using matplotlib

        Does a basic plot of the PMF in matplotlib

        """
        import matplotlib.pyplot as plt
        plt.figure(dpi=300,figsize=(3,3))
        plt.plot(self.center,self.Calc_PMF())
        plt.xlabel("r (angstroms)")
        plt.ylabel("P(r)")
        plt.show()
    



