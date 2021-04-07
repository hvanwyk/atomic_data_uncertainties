import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
"""
NOTE: Currently not used in code
"""


"""
Fit a 'sum of exponentials' function y=h(x,a,b) to a series of 
input-output pairs (xi,yi), where

    h(x,a,b) = sum_{i=1,..,n} a[i]*exp(-b[i]*x) for any x.
    
The problem takes the form of a least squares optimization problem

    min_[a,b] sum_{j=1,...,m} 0.5|h(xi,a,b)-yi|**2
    
The optimization routine is performed by the "curve_fit" method within
the "scipy.optimize" module. To allow for a variable number of terms, 
we need an intermediary 'wrapper function' that transforms the generic
function h into one with a fixed number of parameters.
"""
def wrapper(x, *args):
    """
    Wrapper function in which parameters are passed as a long list. 
    """ 
    n = np.int(len(args)/2)
    a = list(args[0:n])
    b = list(args[n:2*n])
    return h(x, a, b)


def h(x, a, b): 
    """
    Generic model function
    
        h(x,a,b) = sum_{i=0,..,n-1} a[i]*exp(-b[i]*x) for any x.

    Inputs:
    
        x: double, (m,) input vector
        
        a: double, (n,) vector of coefficients
        
        b: double, (n,) vector of scaling parameters
        
    """
    y = np.zeros(len(x))
    for ai, bi in zip(a, b):
        y += ai*np.exp(-bi*x)
    return y


def fit(x,y,a0,b0):
    """
    Determine the optimal parameters a, and b
    """
    n = len(a0)
    p0 = list(a0)
    p0.extend(list(b0))
    popt, pcov = curve_fit(lambda x, *p0: wrapper(x, *p0), x, y, p0=p0, maxfev=100000)
    a = popt[0:n]
    b = popt[n:2*n]
    return a, b
    
def gen_data(x, amplitudes, timeconstants, noise=1): #generate some fake data
    y = np.zeros(len(x))
    for m,t in zip(amplitudes, timeconstants):
        y += m*(1.0-np.exp(-t*x))
    if noise:
        y += np.random.normal(0, noise, size=len(x))
    return y


if __name__=="__main__":
    #
    # Define exact function
    # 
    n = 5       # number of expansion terms
    ae = 3*(np.random.rand(n)-0.5)      # exact coefficients
    be = 0.1*np.random.rand(n)     # exact scalings
    
    #
    # Generate data
    # 
    m = 101  # number of data points
    eps = 1e-2  # noise level
    x = np.linspace(0,10,m)  # x-values    
    y = h(x, ae, be) + eps*np.random.rand(m)  # evaluate exact expansion and pollute
    
    
    for n in np.arange(1,6):
        print(n)
        if n==1:
            a0 = np.ones(n)
            b0 = np.zeros(n)
        else:
            a0 = np.ones(n)
            #a0[0:n-1] = af
            
            b0 = np.zeros(n)
            #b0[0:n-1] = bf
            
        af, bf = fit(x,y,a0,b0)
        plt.plot(x, y, ':', x, h(x,af, bf), '--')
        plt.show()