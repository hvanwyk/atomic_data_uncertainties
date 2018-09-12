import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from mpl_toolkits.mplot3d import Axes3D
import emcee
import corner 

"""
This code snippet illustrates how to wrap a Bayesian sampler around
the mapping from lambda's to the quantities to be used for calibration. 
We want to use the EMCEE package for calibration, which is based on Monte Carlo
sampling. To accommodate a possibly large sample size, we first build a 
surrogate model by creating a linear interpolant over a regular lambda grid.   
In particular, we want to achieve the following.

    1. Interpolate a multivariate mapping f from R^d to R^n, i.e. 
        Mapping x=[x1,x2,...,xd] to y=[f1(x), f2(x), ..., fn(x)]
        based on a regular x-grid. In this implementation, I have split
        the forward mapping into two separate functions (you can combine them
        if you want):
        
            a. fake_autostructure_and_postprocessing works out the quantities
                interest (in this case energies). 
                
            b. fake error function takes the output from the autostructure run
                and computes the deviation from observed values. 
            
    2. Use the python package emcee to generate a MCMC sample of the 
        x vectors in R^d from the posterior distribution, p_post, 
        where 
        
            p_post(x) ~ p_prior(x) p_likelihood(x)
            
    3. Construct a multivariable histogram from the prior
    
    4. Sample from the histogram.
"""

def fake_autostructure_and_postprocessing(x):
    """
    This fake function maps lambda parameters to quantities of interest. 
    
    Input:
    
        x: double, (d, ) array representing a single point in x-space
        
    Output:
    
        y: double, (n, ) array representing a single point in output 
            space.
    """
    #
    # Compute the quantities
    # 
    y1 = np.exp(x[0]**2-2*x[1])*x[1]
    y2 = (x[0]-1)**2 + (x[1]+2)**3
    y3 = np.sin(2*np.pi*x[0])*np.cos(np.pi*x[1])
    
    # 
    # Return output in the form of an (n, ) array
    # 
    return np.array([y1, y2, y3])


def fake_error_function(x, y_nist):
    """
    This fake function computes the error, i.e. the deviation of 
    """
    #
    # Compute the output
    # 
    y = fake_autostructure_and_postprocessing(x)
    
    #
    #  Compute L2 deviation of output from nist   
    # 
    return (y-y_nist)/(1+y_nist)


def log_uniform_prior(x, x_bnd):
    """
    Natural logarithm of the uniform prior, which encodes our initial belief 
    concerning the distribution of x
    
    Inputs: 
    
        x: double, (d,) array of a single point in x space 
        
        x: double, (d,2) array of lower and upper bounds for the x's
        
    
    Output:
    
        log(pdf) = - sum_{i=1,...,d} log(lmd_max_i-lmd_min_i)
    """
    d = len(x)
    x_min, x_max = x_bnd[:,0], x_bnd[:,1]
    
    if np.all((x >= x_min)*(x <= x_max)):
        return -np.sum(np.log(x_max-x_min))
    else:
        return -np.infty
    
    

def log_uniform_likelihood(x, interpolators, y_bnd):
    """
    Logarithm of the likelhood function (using uniform bounds on the error).
    
    Inputs:
    
        x: double, (d,) array of a single point in x space
        
        interpolators: RegularGridInterpolator, n-length list of interpolating
            functions.
            
        y_bnd: double, (d,2) lower and upper bounds on the output y
    """
    n = len(interpolators)
    y_min, y_max = y_bnd[:,0], y_bnd[:,1]
    y = np.zeros(n)
    for i, interpolator in zip(range(n),interpolators):
        y[i] = interpolator(x)
        
    if np.all((y>=y_min)*(y<=y_max)):
        return -np.sum(np.log(y_max-y_min))
    else:
        return -np.infty

    
    
def log_posterior(x, interpolators, x_bnd, y_bnd):
    """
    Compute the log of the posterior density, formed from the logs of the
    likelihood and prior density functions. 
    
    Inputs: 
    
        x: double, (d,) array of a single point in x-space
        
        interpolators: RegularGridInterpolator, n-list of regular grid 
            interpolating functions computed using scipy.RegularGridInterpolator
        
        x_bnd: (d, 2) array of lower and upper bounds for the inputs x 
        
        y_bnd: (n, 2) array of lower and upper bounds for the outputs y
        
    
    Output:
    
        log(post) = log(prior) + log(likelihood) 
        
    NOTE: If you want to be a bit more flexible in the prior and likelihood 
        distributions, you can use flags like 'likelihood_type' or 'prior_type'
        
    """
    lp = log_uniform_prior(x, x_bnd)
    if not np.isfinite(lp):
        return -np.infty
    else:
        return lp + log_uniform_likelihood(x, interpolators, y_bnd)


def interpolators_from_scratch(x_bnd, x_res, n, y_nist):
    """
    Constructs a list of RegularGridInterpolators, one for each output
    
    Inputs:
    
        x_bnd: double, (d,2) array whose ith row contain the lower and upper 
            bounds of the ith input variable.
            
        r_res: int, (d,) array whose ith row contains the number of grid
            points in each dimension.
            
        n: int, dimension of the error vector
        
        y_nist: double, (n, ) array of NIST (or experimental) values 
    """
    #
    # One dimensional grid in each direction
    # 
    d = x_bnd.shape[0]
    X_1D = []
    for i in range(d):
        X_1D.append(np.linspace(x_bnd[i,0],x_bnd[i,1],x_res[i]))
    
    #
    # Multi-dimensional Grid, using meshgrid (ij indexing)
    #
    X = np.meshgrid(*X_1D, indexing='ij')
        
    #
    # Unravel grid into (n_points, d) array
    # 
    x_ravel = np.array([x.ravel() for x in X]).T

    #
    # Generate error data 
    # 
    n_points = x_ravel.shape[0]
    e = np.empty((n_points,n))
    for i in range(n_points):
        x = x_ravel[i,:]
        e[i,:] = fake_error_function(x, y_nist)   
    E = [np.reshape(e[:,j], x_res) for j in range(n)]
    
    #
    # Form the interpolators
    # 
    interpolators = []
    for j in range(n):
        p_fn = RegularGridInterpolator(X_1D, E[j], bounds_error=False)
        interpolators.append(p_fn)
    
    return interpolators


def interpolators_from_data(data_file):
    """
    Constructs the interpolators from data 
    """
    pass


def sample_from_histogram(H, edges, n_samples):
    """
    Generate a random sample from a histogram defined over a hypercube.
    We proceed by conditional sampling, i.e.
    
    1. Pick a random bin based on the histogram values.
    2. Choose a random point inside each subcube. 
    
    """
    # 
    # Generate uniformly distributed random variables
    # 
    U = np.random.rand(n_samples)
    
    # 
    # Normalized cumulative sum from 0 to 1  
    # 
    CDF = np.cumsum(H.ravel())
    CDF = CDF/CDF[-1]
    
    #
    # Select bins randomly.
    # 
    bins = np.searchsorted(CDF, U)
    
    #
    # Determine the bin coordinates within the multi-dimensional grid 
    # 
    n_cells = tuple([len(edge)-1 for edge in edges])
    i_bins = np.array(np.unravel_index(bins, n_cells)).T
    
    d = i_bins.shape[1]
    V = np.random.rand(n_samples, d)
    
    X = np.empty((n_samples,d))
    for j in range(d):
        dxj = np.diff(edges[j])
        i_binsj = i_bins[:,j]
        X[:,j] = edges[j][i_binsj] + V[:,j]*dxj[i_binsj]
    return X

# =============================================================================
# 1. Construct Interpolators
# =============================================================================

# Interval endpoints for each input component
x_bnd = np.array([[0.8,1.2],[0.5,1.3]])

# Resolution in each dimension
x_res = np.array([5,7])

# Fake NIST data
y_nist = fake_autostructure_and_postprocessing(np.array([1,1,1]))

# Construct interpolators
n = len(y_nist) 
interpolators = interpolators_from_scratch(x_bnd, x_res, n, y_nist)


# -----------------------------------------------------------------------------
# Plot interpolators and exact functions
# -----------------------------------------------------------------------------
x0 = np.linspace(x_bnd[0,0], x_bnd[0,1], 31)
x1 = np.linspace(x_bnd[1,0], x_bnd[1,1], 31)
 
X,Y = np.meshgrid(x0,x1)
xy  = np.array([X.ravel(), Y.ravel()]).T

e_exact = []
e_exact.append((np.exp(X**2-2*Y)*Y - y_nist[0])/(1+y_nist[0]))
e_exact.append(((X-1)**2 + (Y+2)**3 - y_nist[1])/(1+y_nist[1]))
e_exact.append((np.sin(2*np.pi*X)*np.cos(np.pi*Y) - y_nist[2])/(1+y_nist[2]))

fig = plt.figure()

#
# Plot Exact surfaces
# 
for j in range(3):
    ax = fig.add_subplot(3,3,j+1, projection='3d')
    ax.plot_surface(X,Y,e_exact[j])
    ax.set_title('Exact Value (Output %d)'%(j+1))
#
# Plot Interpolated Surfaces
# 
for j in range(3):
    ax = fig.add_subplot(3,3,4+j, projection='3d')
    ax.plot_surface(X,Y, interpolators[j](xy).reshape(X.shape))
    ax.set_title('Interpolated Value (Output %d)'%(j+1))
#
# Plot Error
# 
for j in range(3):
    ax = fig.add_subplot(3,3,7+j, projection='3d')
    ax.plot_surface(X,Y, interpolators[j](xy).reshape(X.shape)-e_exact[j])
    ax.set_title('Interpolation Error (Output %d)'%(j+1))
    
plt.show()


# =============================================================================
# Run MCMC Code
# =============================================================================
#
# Likelihood is based on 5% error in each component 
# 
y_bnd = 0.05*np.array([[-1,1],[-1,1],[-1,1]])

#
# Specify the dimension of the input space and the number of starting points
#
n_dim, n_walkers = 2, 100

#
# Specify starting points for each Markov chain (in a tight ball around optimum)
#
pos = [np.array([1,1]) + 1e-4*np.random.randn(n_dim) for i in range(n_walkers)]

#
# Initialize the sampler
#
sampler = emcee.EnsembleSampler(n_walkers, n_dim, log_posterior, 
                                args=(interpolators, x_bnd, y_bnd))
#
# Run the MCMC routine
#
sampler.run_mcmc(pos, 1000);

#
# The sampler.chain has shape (n_walkers, n_steps, n_dim) = (100, 1000, 2)
# Reshape into a (100'000,2) array of samples (but throw away initial sample -> (95'000,2)
#
samples = sampler.chain[:, 50:, :].reshape((-1, n_dim))

# -----------------------------------------------------------------------------
# Visualize the posterior density
# -----------------------------------------------------------------------------
corner.corner(samples, labels=["$\lambda_1$", "$\lambda_2$"], truths=[1,1])
plt.show()

# -----------------------------------------------------------------------------
# Plug the Monte Carlo samples from the posterior into the function and comput
# a histogram for the output. Make sure the point [0,0,0] has nonzero density   
# -----------------------------------------------------------------------------
n_samples = samples.shape[0]
E_samples = np.empty((n_samples,n))
for i in range(n_samples):
    E_samples[i,:] = fake_error_function(samples[i,:], y_nist)

corner.corner(E_samples, labels=["$E_1$","$E_2$","$E_3$"], truths=[0,0,0])
plt.show()


# =============================================================================
# Compute histogram
# =============================================================================
H, edges = np.histogramdd(samples, bins=(20,20), normed=True)
CH = np.cumsum(H.ravel())
plt.plot(CH)
plt.show()


# =============================================================================
# Sample from histogram
# =============================================================================
y = sample_from_histogram(H, edges, 10000)

corner.corner(y, labels=["$\lambda_1$", "$\lambda_2$"], truths=[1,1] )
plt.show()
