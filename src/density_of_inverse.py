# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 13:16:13 2017

@author: hans-werner

NOTE: We no longer use the inverse density estimation in our work.
"""

# ============================================================================
# Imports 
# ============================================================================
import numpy as np
#from numpy.linalg import norm
import matplotlib.pyplot as plt
import math
from itertools import chain, combinations, product
from numbers import Number 
import scipy.sparse as sp
from mpl_toolkits.mplot3d import Axes3D # @UnresolvedImport

# ============================================================================
# Classes
# ============================================================================
class BlockFunction:
    """
    Function defined on a single Rectangular block
    
    Attributes:
    
     
    
    Methods:
    
    """
    def __init__(self, x_bnd, f, dfdx, address=None):
        """
        Constructor 
        
        Inputs:
        
            x_bnd: double, (n,2) array of coordinate endpoints specifying block
            
            f: function, to be approximated OR
               double, function value at block midpoint
            
            dfdx: function, gradient of f OR
                  (n,) array of directional derivatives at the block midpoint
            
            address [None]: int, address of Block (for use in GridFunction).
        """
        if len(x_bnd.shape) == 1:
            #
            # One dimensional case
            # 
            x_mid = 0.5*(x_bnd[0]+x_bnd[1])
            self.__is_1d = True
        else:
            #
            # Multidimensional case
            #
            x_mid = 0.5*(x_bnd[:,0]+x_bnd[:,1])
            self.__is_1d = False
        self.__x_mid = x_mid
        n = x_bnd.shape[0]
        #
        # Function may be explicit or a number
        # 
        if callable(f):
            self.__f0 = f(self.__x_mid)
        elif isinstance(f, Number):
            self.__f0 = f
        else:
            raise Exception('f is a function or a number.') 
        #
        # Gradient may be a function or a vector
        # 
        if callable(dfdx):
            #
            # gradient given as function, evaluate it
            # 
            self.__dfdx = dfdx(self.__x_mid)
        elif type(dfdx) is np.ndarray:
            #
            # gradient value specified, use it
            # 
            assert len(dfdx)==n, 'Gradient dimension mismatch.'
            self.__dfdx = dfdx
        elif dfdx is None:
            #
            # gradient not specified, approximate it by average gradient
            #
            dfdx = np.empty(n)
            I = np.eye(n)
            x0 = x_bnd[:,0]
            for i in range(n):
                dxi = x_bnd[i,1]-x_bnd[i,0]
                dx = I[:,i]*dxi  # change in x[i]
                df = f(x0+dx)-f(x0)  # change in f
                dfdx[i] = df[0]/dxi  # average rate of change
            self.__dfdx = dfdx
        else:
            raise Exception('dfdx is a function or a vector')
            
        self.__n = x_bnd.shape[0]
        self.__x_bnd = x_bnd
        self.__address = address
        self.__p = 0
   
    
    def bnd(self):
        """
        Return boundary delimiter
        """
        return self.__x_bnd
           

    def dfdx(self):
        """
        Return gradient
        """
        return self.__dfdx
 
 
    def f(self):
        """
        Return function value at midpoint
        """
        return self.__f0        
        
        
    def midpoint(self):
        """
        Return box midpoint
        """        
        return self.__x_mid
        
        
    def address(self):
        """
        Return block's address
        """
        return self.__address
        
        
    def assign_address(self, address):
        """
        Assign an identifying address to the block        
        """
        self.__address = address
        
    
    def probability(self):
        """
        Return the probability of the cell
        """
        return self.__p        
        
        
    def assign_probability(self, p):
        """
        Assign a probability to a given cell
        """
        self.__p = p        
        
        
    def volume(self):
        """
        Compute the volume of a cell
        """
        xb = self.bnd()
        if self.__is_1d:
            return xb[1]-xb[0]
        else:
            return np.prod(xb[:,1]-xb[:,0])
                
        
    def extrema(self, min_or_max='both'):
        """
        Compute the local extrema of the piecewise linear function over the
        hyperrectangle.
        
        Outputs: 
        
            x_min: double, (n,) argmin of f over block
            
            f_min: double, min of f over block
        """
        xv = self.vertices()  # extract block vertices        
        Lv = self.evaluate(xv)  # evaluate tangent plane at vertices
        if min_or_max == 'min':
            i_min = np.argmin(Lv)
            x_min = xv[i_min] if self.__is_1d else xv[:,i_min]
            return x_min, Lv[i_min]
        elif min_or_max == 'max':
            i_max = np.argmin(Lv)  # get index of first occurence of maximum
            x_max = xv[i_max] if self.__is_1d else xv[:,i_max]
            return x_max, Lv[i_max]
        elif min_or_max == 'both': 
            i_min = np.argmin(Lv)
            i_max = np.argmax(Lv)
            x_min = xv[i_min] if self.__is_1d else xv[:,i_min]
            x_max = xv[i_max] if self.__is_1d else xv[:,i_max]
            return x_min, Lv[i_min], x_max, Lv[i_max]
        elif min_or_max == 'range':
            return Lv[i_min], Lv[i_max]
        else: 
            raise Exception('Variable min_or_max can be "min", ' + \
                            '"max", "range", or "both".')
               
        
    def evaluate(self, x):
        """
        Evaluate the piecewise linear approximation of f at a collection of 
        points x = [x1, ..., xk]
        
        Output:
        
            L: double (k, )  f0 + dfdx*(x-x0)
        """
        if self.__is_1d:
            L = self.__f0 + self.__dfdx*(x-self.__x_mid)
            #
            # Assign function value 'nan' to points outside interval
            # 
            in_box = self.in_box(x)
            if isinstance(x, Number):
                L = L if in_box else np.nan
            else:
                L[~in_box] = np.nan
        else:
            #
            # Check if dimensions match
            # 
            assert x.shape[0] == self.__n, \
                'Dimension of point differs from that of box.'
            
            # 
            # Turn into a column vector if necessary
            # 
            if len(x.shape) == 1:
                x = np.reshape(x, (-1,1))    
            k = x.shape[1]
        
            #
            # Evaluate hyperplane
            # 
            dx = x - np.tile(self.__x_mid, (k,1)).transpose()
            L = self.__f0 + np.dot(self.__dfdx, dx)
        
            # 
            # Assign nan to function values of points outside box
            #
            in_box = self.in_box(x)
            L[~in_box] = np.nan

        return L
             
                
    def in_box(self, x, tol=1e-10): 
        """
        Determine which points in x = [x1, x2, ..., xk] lie in the block
        
        Input: 
        
            x: double, (n,k) array of points 
            
        
        Output:
        
            in_box: bool, (k,) vector with components True or False, depending
                on whether they lie inside the box
        """
        bnd = self.bnd()
        if self.__is_1d:
            #
            # One dimensional case
            #
            if isinstance(x, Number):
                in_box = (x >= bnd[0]-tol) and (x<=bnd[1]+tol)
            else:
                assert len(x.shape) == 1, 'Input should be a number/vector.'
                in_box = (x >= bnd[0]-tol)*(x<=bnd[1]+tol)
        else:
            # 
            # Multidimensional case
            #
            assert x.shape[0] == self.__n, \
                'Dimension of point differs from that of box.'
            
            if len(x.shape) == 1:
                k = 1
            else:
                k = x.shape[1]
            
            in_box = np.empty(k, dtype=bool)
            for j in range(k):
                in_box[j] = np.all( ( x[:,j] >= bnd[:,0]-tol ) * \
                                    ( x[:,j] <= bnd[:,1]+tol )  )
        return in_box
    

    
    def vertices(self):
        """
        Return block vertices 
        
        Output:
        
            xv: double, (n,2^n) array of vertices
        """
        if self.__is_1d: 
            return self.__x_bnd
        else:
            xv = np.array([p for p in product(*tuple(self.__x_bnd))])
            return xv.transpose()
    

    def slab_vol(self, f_min, f_max):
        """
        Compute the volume of a slab of contour lines within a block, defined
        by a range of function values.
        
        Inputs:
        
            f_min, f_max: double, minimum and maximum function values
            
        Output:
        
            slab_volume: double, volume of slab within block.
        """
        #
        # Hyperplane parameters
        #
        bnd = self.bnd()
        dfdx = self.dfdx()
        f0 = self.f()
        x0 = self.midpoint()
      
        if self.__is_1d:
            #
            # One dimensional case
            #
            if abs(dfdx) < 1e-12:
                # Horizontal 
                vol = (bnd[1]-bnd[0]) if (f_min <= f0 and f0<=f_max) else 0
            else:
                x_inv = [(f_min-f0)/dfdx+x0, (f_max-f0)/dfdx+x0]
                x_min, x_max = tuple(sorted(x_inv))
                #print('Going Back')
                #print(self.evaluate(x_min))
                #print(self.evaluate(x_max))
                #print('f^(-1)[%.3f,%.3f]=[%.3f,%.3f]'%(f_min,f_max,x_min,x_max))
                vol = np.max([np.min([x_max-bnd[0], bnd[1]-x_min]), 0]) 
        else:
            #
            # Exctract degenerate directions
            #
            dx = bnd[:,1] - bnd[:,0]
            i_dgn = abs(dfdx) < 1e-12
            dfdx = dfdx[~i_dgn]
            x0 = x0[~i_dgn]
            z_min = f_min - f0 + dfdx.dot(x0-bnd[~i_dgn,0]).squeeze()
            z_max = f_max - f0 + dfdx.dot(x0-bnd[~i_dgn,0]).squeeze()
            v = self.vol_halfspace_unitcube(dfdx*dx[~i_dgn], z_max) - \
                self.vol_halfspace_unitcube(dfdx*dx[~i_dgn], z_min )
            vol = np.prod(dx)*v
            
        return vol
        
        
    def vol_halfspace_unitcube(self, w, z):
        """
        Compute the volume of the intersection of a halfspace <w,x> <= z with 
        the unit hypercube [0,1]^n in R^n.
        
        Inputs:
            
            z: double, z intercept of hyperplane
            
            w: double, (n,) normal vector of generating hyperplane with no non-zero
                entries.
            
        
        Outputs:
        
            V: double >0, volume 
            
            
        Reference: 
        
            "SLICES, SLABS, AND SECTIONS OF THE UNIT HYPERCUBE" by 
            JEAN-LUC MARICHAL AND MICHAEL J. MOSSINGHOFF.
        
        """
        if any(np.abs(w) < 1e-12):
            #
            # Degenerate hyperplane
            #
            return 
        #
        # Check for trivial intersections
        #
        n = len(w)
        bnd = np.array([np.zeros(n),np.ones(n)]).transpose()
        vertices = np.array([p for p in product(*tuple(bnd))]).transpose()
        z_test = w.dot(vertices)
        if all(z_test > z):
            #
            # Halfspace does not intersect with hypercube
            #
            return 0
        elif all(z_test <= z):
            #
            # Hypercube fully contained in halfspace
            # 
            return 1
            
        # 
        # Non-trivial intersection
        #
        vol = 0
        n = len(w)
        I = set(range(n))
        for K in chain.from_iterable(combinations(I,r)\
            for r in range(len(I)+1)):
            # 
            # Iterate over power set of [1,...,n]
            #
            vol += (-1)**len(K)*np.max([z-np.sum(w[sorted(K)]),0])**n            
        vol = vol/math.factorial(n)/np.prod(w)
        return vol


class GridFunction():
    """
    Description: A piecewise approximation of a function f, obtained by 
        imposing a grid on the domain and linearizing f at the midpoint of 
        each block.
    
    
    Attributes: 
    
        __grid: BlockFunction, list of hyperplanes that approximate the 
            function f within each sub-block. 
            * keys: tuples specifying address of block in grid
            * values: BlockFunction object
            
        __n: int, dimension
        
        __n_intervals: int, n-tuple the number of subdivisions in each 
            direction.
            
        __bnd: double, (n,2) endpoints of the hypercube domain of GridFunction.
    """
    def __init__(self, bnd, n_intervals, f, dfdx=None):
        """
        Constructor
        
        Inputs:
        
            bnd: double, (n,2) array containing the hypercube boundary.        
        
            n_intervals: int, (n,) tuple containing the number of subdivisions 
                in each direction.
                
            f: function, to be approximated OR vector of (n_cells,) function 
                values at the midpoint of each cell.
            
            dfdx: function, (n,) valued gradient function OR (n_cells,n) array
                of gradients at the midpoint of each cell.
        """
        
        if isinstance(n_intervals, Number):
            #
            # One dimensional
            #
            n = 1
            dx = (bnd[1]-bnd[0])/n_intervals
            grid = {}
            for address in range(n_intervals):
                block_bnd = np.array([bnd[0]+address*dx, \
                                      bnd[0]+(address+1)*dx])
                grid[address] = BlockFunction(block_bnd, f, dfdx, address)        
        else:
            #
            # Multidimensional case
            #
            assert type(n_intervals) is tuple,\
                'Store the number of subdivisions in a tuple.'
            n = len(n_intervals)
            dx = (bnd[:,1]-bnd[:,0])/np.array(n_intervals)
            grid = {}
            range_tuple = tuple([range(ni) for ni in n_intervals])
            for address in product(*range_tuple):
                block_bnd = np.array( [bnd[:,0] + np.array(address)*dx, \
                                       bnd[:,0] +(np.array(address)+1)*dx] )
                block_bnd = block_bnd.transpose()                                
                grid[address] = \
                    BlockFunction(block_bnd, f, dfdx=dfdx, address=address)       
        
        # 
        # Index blocks numerically
        # 
        i2a = {}
        i = 0
        for address in grid.keys():
            i2a[i] = address
            i += 1
        
                
        self.__grid = grid
        self.__n_cells = i
        self.__i2a = i2a
        self.__n = n
        self.__n_intervals = n_intervals
        self.__bnd = bnd
        self.__has_histogram = False
        self.__min_tuple = None
        self.__max_tuple = None
        self.__p_range = None
        self.__output_histogram = None
    
    def block(self, i):
        """
        Return the ith block
        """
        return self.__grid[self.__i2a[i]]
        
        
    def grid(self):
        """
        Return the grid
        """             
        return self.__grid
        
        
    def dim(self):
        """
        Return the dimension of the problem
        """
        return self.__n
        
        
    def n_intervals(self):
        """
        Return a tuple of the number of subdivisions in each direction
        """
        return self.__n_intervals
        
        
    def bnd(self):
        """
        Return the left and right endpoints of the hypercube 
        """
        return self.__bnd


    def enclosing_block(self, x):
        """
        Find a block that contains the given (set of) point(s) x, None otherwise
        
        Inputs:
                
            x: double, (n, k) array of points for which we are searching the
                enclosing block.
        
        TODO: 
        """
        pass
    
    
    def global_extrema(self, min_or_max='both'):
        """
        Compute the maximum/minimum of the piecewise linear approximation of 
        q over an n-dimensional rectangular grid. 
        
        Input:
            
            min_or_max: str, 'min', 'max', 'both', or 'range' 
                Note: If 'min', 'max', or 'both', the function returns the 
                    location, function value, and block containing the global
                    extrem(um/a). If 'range', the function only returns the
                    minimum and maximum function values.
                
        
        Output:
        
            x_min, f_min, b_min: double, location, function value and block 
                that contains the global min. (min_or_max='min' or 'both')
            
            x_max, f_max, b_max: double, locations, function value and block
                that contains the global max. (min_or_max='max' or 'both')
        
        Note: The GridFunction may attains its extrema at multiple locations.
            This method only returns one such location.
        """
        if self.__min_tuple is not None:
            assert self.__max_tuple is not None, \
                'The global min and max should have been computed together.'
            x_min, f_min, b_min = self.__min_tuple
            x_max, f_max, b_max = self.__max_tuple
            
        else:
            x_min, f_min, b_min = None, None, None
            x_max, f_max, b_max = None, None, None
            #
            # Iterate over blocks
            #
            for block in self.grid().values():
            
                # Compute local extrema
                x_min_loc, f_min_loc, x_max_loc, f_max_loc = \
                    block.extrema(min_or_max='both')
                
                # Update global minimum if necessary
                if (f_min is None) or f_min_loc < f_min:
                    x_min, f_min, b_min  = x_min_loc, f_min_loc, block
                
                if (f_max is None) or f_max_loc > f_max:
                    x_max, f_max, b_max = x_max_loc, f_max_loc, block
                    
        self.__min_tuple = (x_min, f_min, b_min)
        self.__max_tuple = (x_max, f_max, b_max)
        #
        # Decide what to return
        #
        if min_or_max == 'min':
            return x_min, f_min, b_min
        elif min_or_max == 'max':
            return x_max, f_max, b_max
        elif min_or_max == 'both':
            return x_min, f_min, b_min, x_max, f_max, b_max
        elif min_or_max == 'range':
            return f_min, f_max
        else:
            raise Exception('Use either "min", "max", "range", or ' + \
                            '"both" for input "min_or_max".')


    def monotone_path(self, x, address, direction='up'):
        """
        Compute the monotone path transverse to the function's contours along 
        the path of steepest (i) increase (gradient), or (ii) decrease 
        (-gradient). 
            
        Inputs:
                
            x: x-value at which path begins
            
            address: int, n-tuple address of block in which initial point lies.
            
            direction: str, specifying whether the path goes along the gradient
                ('up') or against it ('down')
     
        
        Outputs: 
        
            path_x: double, list of (n,2) arrays containing the initial and final 
                x-values within each block. 
            
            path_f: double, list of (2,) arrays containing the approximate function 
                values at the points in 'path_x'
            
        """ 
        block = self.grid()[address]
        path_x = []
        path_f = []
        while True:
            fx = block.evaluate(x)
            block_bnd = block.bnd()
            if direction == 'down':
                p = -block.dfdx()
            elif direction == 'up':
                p = block.dfdx()            
            else:
                raise Exception('Use either "up" or "down" for "direction".')
            
            #
            # Destination boundary (to which p points)
            # 
            dest_bnd = 0.5*(block_bnd[:,0]+block_bnd[:,1]) + \
                       0.5*np.sign(p)*(block_bnd[:,1]-block_bnd[:,0])

           
            #
            # Compute stepsize to nearest boundary edge
            # 
            dx = dest_bnd - x
            a_vec = np.where(p!=0, np.abs(dx/p), np.inf)
            i = np.argmin(a_vec)
            a = a_vec[i]
            
            if np.abs(a) < 1e-13:
                #
                # Local min/max reached -> done!
                # 
                return path_x, path_f
        
            #
            # Update x, fx -> store pairs in path, f
            # 
            x_new = x + a*p
            fx_new = block.evaluate(x_new)
            path_x.append(np.array([x,x_new]).transpose())
            path_f.append(np.array([fx,fx_new]))
            
            #
            # Move to next block
            #
            address = np.array(block.address())
            address[i] += int(np.sign(p[i]))
            if any(address < 0):
                #
                # Reached boundary
                #
                return path_x, path_f
            
            block = self.grid()[tuple(address)]
            x = x_new
            fx = fx_new        


    def set_output_histogram(self, f_grid, f_prob):
        """
        Define the histogram for the output.
        
        Inputs: 
        
            f_grid: double, (k+1, ) ordered vector of f-values between f_min
                and f_max.
                
            f_prob: double >0, (k,) vector of probabilities associated with
                each subdivision.
        """
        
        #
        # Checks
        # 
        assert np.abs(f_prob.sum()-1) < 1e-12, \
            'Probabilities should add up to 1.'
            
        assert len(f_grid) == len(f_prob)+1, \
            'Number of subintervals and length of histogram incompatible.'
        #
        # Assignment
        #            
        self.__output_histogram = (f_grid, f_prob)        


    def output_histogram(self):
        """
        Return the grid and probabilities of the histogram of the output.
        """
        return self.__output_histogram

        
    def compute_histogram(self):
        """
        Approximate the density function of the input space, by means of the 
        inverse sensitivity method. 
        """
        if self.output_histogram() is None:
            raise Exception('Specify histogram for output.')
        else:
            f_grid, f_prob = self.output_histogram()
        #
        # Ensure the range of f_grid lies within f_min and f_max
        #         
        f_min, f_max= self.global_extrema(min_or_max='range')
        if f_grid.min() > f_min:
            f_grid = np.hstack((f_min, f_grid))
            f_prob = np.hstack((0, f_prob))
            
        if f_grid.max() < f_max:
            f_grid = np.hstack((f_grid, f_max))
            f_prob = np.hstack((f_prob, 0))
        
        n_partitions = len(f_prob)
        n_cells = self.__n_cells
        #
        # Initialize sparse matrix of conditional probabilities
        #            
        rows = np.empty(n_cells*n_partitions)
        cols = np.empty(n_cells*n_partitions)
        vals = np.empty(n_cells*n_partitions)
        count = 0
        
        volA = np.zeros(n_partitions)
        for j in range(n_partitions):                        
            for i in range(n_cells):                 
                block = self.block(i)
                
                # Intersection of block with inverse image of f-range
                block_finvA = block.slab_vol(f_grid[j],f_grid[j+1])  
                                
                if np.abs(block_finvA) > 1e-16:
                    rows[count] = i
                    cols[count] = j
                    vals[count] = block_finvA
                    count += 1
                
                # Update total area
                volA[j] += block_finvA
                
        
        rows = rows[:count]
        cols = cols[:count]
        vals = vals[:count]
        
        V = sp.coo_matrix((vals,(rows,cols)),shape=(n_cells,n_partitions), \
                          dtype=np.float)
        #print(V.toarray())
        #print(volA)
        #print(np.array(V.sum(axis=0)).ravel()-volA)
                          
        #
        # Row normalize
        # 
        V = V.transpose().tocsr()
        row_sums = np.array(V.sum(axis=1)).ravel()
        ri, _ = V.nonzero()
        V.data = V.data/row_sums[ri]
        
        #print(V.toarray())

        #
        # Compute cellwise probabilities
        #
        P = V.transpose().dot(f_prob)
        PP = np.empty(self.__n_intervals)
        #print(P.sum())
        for i in range(n_cells):
            block = self.block(i)
            block.assign_probability(P[i])
            PP[block.address()] = P[i]
        self.__has_histogram = True
        self.__P = PP
        
   
    def probability_range(self):
        """
        Return and store the minimum and maximum cellwise probability values 
        (for plotting).
        """
        if self.__p_range is None:
            if not self.__has_histogram:
                raise Exception('Histogram has not been computed.')
            
            # Initialize extrema
            p_min, p_max = 1, 0
            for block in self.grid().values():
                p_block = block.probability()
                if p_block < p_min:
                    p_min = p_block
                    
                if p_block > p_max:
                    p_max = p_block
                
            self.__p_range = (p_min, p_max)
        return self.__p_range
        

    def plot(self, what):
        """
        Plots
        
            what: str, what to plot ('function', 'histogram', )
        """
        assert self.dim() <= 2, 'Can only visualize in 1D or 2D.'
        fig = plt.figure()
        if what == 'function':
            if self.dim()==1:
                ax = fig.gca()
                for block in self.grid().values():
                    x = block.vertices()
                    y = block.evaluate(x)
                    ax.plot(x,y,'-k')
            else:   
                f_min, f_max = self.global_extrema(min_or_max='range')
                ax = fig.gca(projection='3d')
                for block in self.grid().values():
                    x = block.vertices()
                    [X,Y] = np.meshgrid(x[0,:],x[1,:])
                    xy_inputs = np.array([X.ravel(),Y.ravel()])
                    Z = block.evaluate(xy_inputs).reshape(X.shape)
                    ax.plot_surface(X,Y,Z, cmap='viridis', \
                                    vmin=f_min, vmax=f_max)
        elif what == 'histogram':
            if self.dim()==1:
                ax = fig.gca()
                for block in self.grid().values():
                    x = block.vertices()
                    p = block.probability()*np.ones(x.shape)
                    ax.fill_between(x,0,p,facecolor='blue', edgecolor='k', \
                                    alpha=0.5)
            else:
                ax = fig.gca(projection='3d')
                p_min, p_max = self.probability_range()
                for block in self.grid().values():
                    x = block.vertices()
                    [X,Y] = np.meshgrid(x[0,:],x[1,:])
                    Z = np.ones(X.shape)*block.probability()
                    ax.plot_surface(X,Y,Z,alpha=0.5, cmap='viridis', \
                                    vmin=p_min, vmax=p_max)
    
        elif what == 'contour':
            if self.dim() ==1:
                pass
            else:
                ax = fig.gca()
                f_min, f_max = self.global_extrema(min_or_max='range')
                f_grid, f_prob = self.output_histogram()
                for block in self.grid().values():
                    bnd = block.bnd()
                    x = np.linspace(bnd[0,0],bnd[0,1],20)
                    y = np.linspace(bnd[1,0],bnd[1,1],20)
                    [X,Y] = np.meshgrid(x,y)
                    xy = np.array([X.ravel(),Y.ravel()])
                    Z = block.evaluate(xy).reshape(X.shape)
                    c = ax.contourf(X,Y,Z,f_grid, cmap='viridis', origin='lower',\
                                vmin=f_min, vmax=f_max, alpha=0.5)
                plt.colorbar(c)
        plt.show()
                    