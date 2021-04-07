import numpy as np
from scipy.stats import rv_continuous
from scipy.stats import multivariate_normal
from scipy.stats import norm
from scipy.stats import truncnorm
from scipy.stats import uniform
from scipy.stats import beta

#
# TODO: Move to bayesian_methods
# 
def log_likelihood(x, y_obs, interpolators, rvs, error_type='mixed'):
    """
    We assume the measurements were made independently, so that
    
    
    Inputs:
    
        x: double, (d,) array representing a single point in input 
            space.
            
        y_obs: double, (n,) array representing the measured observations
            of the outputs. 
            
        interpolators: RegularGridInterpolator, length n list of interpolating
            functions used to approximate the forward map f:x -> y
            
        rvs: rv_continuous, length n list of random variables with specified 
            densities. If rvs is a single instance of an "rv_continuous" object,
            we assume all n output random variables have the same distribution.
        
        error_type: s(yi-yi_obs)/np.abs(yi_obs)tr, type of univariate error estimate 
            'absolute' (yi-yi_obs), 
            'relative' (yi-yi_obs)/|yi_obs|, or 
            'mixed'    (yi-yi_obs)/(1+|yi_obs|). 
        
    Outputs:
    
        log_pdf = sum_{i=1,...,n} log_pi[ fi(x)-yi_obs ]
    """
    n = len(y_obs)
    
    #
    # Check that there are enough interpolators
    #
    assert len(interpolators)==n, 'The number of interpolators should equal '+\
        'the dimension of observation space.'
    
    if type(rvs) is list:
        # 
        # Possibly different random variables for each output
        # 
        # Check that there are n of them. 
        assert len(rvs)==n, 'Number of random variable incompatible with '+\
            'dimension of the output/observation space.'
    else:
        #
        # Single random variable -> make n copies
        # 
        assert isinstance(rvs, rv_continuous), 'Input "rvs" should be a '+\
        '(list of type) "rv_continuous".'
        rvs = [rvs]*n
    #
    # Compute log-pdf of joint distribution
    # 
    log_pdf = 0
    for i in range(n):
        #
        # Compute the ith model output
        #
        yi = interpolators[i](x)
        
        #
        # Get ith observation 
        #
        yi_obs = y_obs[i]
        
        #
        # Compute the error
        #
        if error_type=='absolute':
            ei = yi-yi_obs
        elif error_type=='relative':
            ei = (yi-yi_obs)/np.abs(yi_obs)
        elif error_type=='mixed':
            ei = (yi-yi_obs)/(1+np.abs(yi_obs))
        
        #
        # Compute the log-pdf of the ith random variable
        # 
        log_pdf += rvs[i].logpdf(ei)
        
    return log_pdf


def log_prior(x, rvs):
    """
    Compute the log prior density on the input parameters
    
    Inputs:
    
        x: (d,) array representing a single point in input space
        
        density_type: str, name of density, 'mvn',
        
        density_parameters: dict, of 
    """
    d = len(x)
    #
    # Check random variables
    # 
    if type(rvs) is list:
        assert len(list)==d, 'Number of random variables not compatible '+\
            'with size of input vector.'
        
        for rv in rvs:
            assert isinstance(rv, rv_continuous), 'At least one entry in '+\
                'input "rvs" is not of type "rv_continuous".'
    else:
        assert isinstance(rvs, rv_continuous), 'Input "rvs" should be an '+\
            'instance of object "rv_continuous".'
        if isinstance(rvs, multivariate_normal):
            #
            # Special treatment for multivariate normal distribution 
            # 
            pass
        else:
            #
            # List 
            #
            rvs = [rvs]*d
        
    #
    # Compute log prior
    # 
    if isinstance(rvs, multivariate_normal):
        #
        # Multivariate normal density (with possible correlations)
        # 
        return rvs.logpdf(x)
    else:
        #
        # Multivariate density of independent random variables.
        # 
        log_prior = 0
        for i in range(d):
            #
            # ith component of vector
            # 
            xi = x[i]
            
            #
            # ith random variable
            # 
            rvi = rvs[i]
            
            #
            # Compute log prior
            # 
            log_prior += rvi.logpdf(xi)
        return log_prior    
    
    
def log_posterior(x, y_obs, interpolators, prior_rvs, 
                  likelihood_rvs, error_type):
    """
    Compute the log posterior of x, given y_obs and a statistical model 
    for the error.
    """
    lp = log_prior(x, prior_rvs)
    if not np.isfinite(lp):
        #
        # Log prior is -oo, no need to compute likelihood
        # 
        return -np.infty
    else:
        #
        # Add log likelihood to log prior
        # 
        ll = log_likelihood(x, y_obs, interpolators, \
                            likelihood_rvs, error_type=error_type)
        return ll + lp