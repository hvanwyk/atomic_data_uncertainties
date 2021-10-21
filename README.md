

# Uncertainties in Atomic Data (uad)
Systematic treatment of uncertainties in atomic structure.The first development code is for Dielectronic recombination (DR) and is in the recombination directory. 


For recombination, the method uses a Bayesian analysis Markov-Chain Monte-Carlo approach to determine a distribution function of orbital scaling parameters. This distribtion function is then carried through a DR calculation to produce uncertainties on total DR rate coefficients.


## Structure of the code


## Software Dependencies

- ```emcee```: Bayesian analysis package ```emcee``` to estimate the distribution of the input parameters.
- ```ColRadPy```: Collisional radiative models
- ```pyatomdb```: Python interface for downloading astrophysical data from the ATOMDB
- ```Autostructure```: 
