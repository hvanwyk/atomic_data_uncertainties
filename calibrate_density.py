# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 13:16:13 2017

@author: hans-werner

"""

# ============================================================================
# Imports 
# ============================================================================
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
import math
from itertools import chain, combinations, product
from numbers import Number 
import scipy.sparse as sp
from mpl_toolkits.mplot3d import Axes3D

# ============================================================================
# Classes
# ============================================================================
class BlockFunction:
    """
    Function defined on a single Rectangular block
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
        else: 
            raise Exception('Variable min_or_max can be "min", ' + \
                            '"max", or "both".')
               
        
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
            
            min_or_max: str, 'min', 'max' or 'both' if min/max/or both extrema
                are sought. Note, all quantities are computed regardless of 
                this string. 
                
        
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
        else:
            raise Exception('Use either "min", "max", or "both" for ' + \
                            'input "min_or_max".')


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
        x_min, f_min, b_min, x_max, f_max, b_max = self.global_extrema()
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
        ri, ci = V.nonzero()
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
                x_min, f_min, b_min, x_max, f_max, b_max = \
                    self.global_extrema()
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
                x_min,f_min,b_min,x_max,f_max,b_max = self.global_extrema()
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
      
                


# ============================================================================
# Tests
# ============================================================================
def test_01():
    """
    General Test for GridFunction over a 2D square
    """
    # Function and gradient 
    f = lambda x: x[0]*x[1]*np.exp(-(x[0]**2+1.25*x[1]**2-1))
    df_dx = lambda x: \
        np.array([ x[1]*np.exp(-(x[0]**2+1.25*x[1]**2-1))*(1-2*x[0]**2), \
                   x[0]*np.exp(-(x[0]**2+1.25*x[1]**2-1))*(1-2.5*x[1]**2)])
    # Hypercube range               
    bnd = np.array([[0,2],[0,2]]) 
    
    # Refinement level
    n_intervals = (15,15)
    
    # Construct gridfunction
    g = GridFunction(bnd, n_intervals, f, df_dx)
    
    # Subdivide range space into partitions
    x_min, f_min, b_min = g.global_extrema('min')
    x_max, f_max, b_max = g.global_extrema('max')
    n_partitions = 10
    f_grid = np.linspace(f_min, f_max, n_partitions+1)
    
    # Impose uniform density
    f_prob = np.ones(n_partitions)/n_partitions
    
    # Compute histogram on grid
    g.set_output_histogram(f_grid, f_prob)
    g.compute_histogram() 
    """
    n_cells = np.prod(n_intervals)
    p = np.empty(n_cells)
    count = 0
    for block in g.grid().values():
        p[count] = block.probability()
        count += 1
    assert abs(p.sum()-1) < 1e-12, 'Probabilities should add up to 1.'
    """
    # Plot results
    g.plot('function')
    g.plot('histogram')
    g.plot('contour')
    
def test_02():
    """
    Test the computation of volumes of slabs
    FIXME: vol_halfspace_unitcube is now a BlockFunction method.
    """    
    # ------------------------------------------------------------------------
    # Test: vol_halfspace_unitcube
    # ------------------------------------------------------------------------
    block = BlockFunction(np.array([0,1]), lambda x: 1, lambda x: 0)
    
    # 2D
    w = np.array([-2,-1])
    z = -2
    v = block.vol_halfspace_unitcube(w,z)
    assert abs(v-0.25) < 1e-12, 'Volume should be 1/4.'
    
    v = block.vol_halfspace_unitcube(-w,-z)
    assert abs(v-0.75) < 1e-12, 'Volumes should add up to 1'
    
    # 3D
    w = np.array([0,1,0])
    z = 1
    v = block.vol_halfspace_unitcube(w,z)
    assert v is None, 'Degeneracy, answer should be None.'
    
    w = np.array([1,1,1])
    z = 1
    v = block.vol_halfspace_unitcube(w,z)
    assert abs(v-1/6) < 1e-12, 'Volume should be 1/6.'
    
    # ------------------------------------------------------------------------
    # Test slab_vol
    # ------------------------------------------------------------------------
    
    # Horizontal hyperplane: unit hypercube
    f1 = lambda x: x[1]
    bnd = np.array([[0,0,0],[1,1,1]]).transpose()
    df1_dx = lambda x: np.array([0,1,0]).transpose()
    bf = BlockFunction(bnd, f1, df1_dx)
    v = bf.slab_vol(0.5,1)
    assert abs(v-0.5)< 1e-12, 'Volume should be 1/2.'
    
    # Horizontal hyperplane, nonstandard hypercube
    bnd = np.array([[0,1,0],[0.5,2,2]]).transpose()
    bf = BlockFunction(bnd, f1, df1_dx)
    assert abs(bf.slab_vol(-1,1)) < 1e-12, 'Volume should be 0'
    assert abs(bf.slab_vol(1,4)-1) < 1e-12, 'Volume should be 1'     
    
    # Skew hyperplane
    f2 = lambda x: x[0] + x[1] - 2
    df2_dx = lambda x: np.array([1,1])
    bnd = np.array([[1,1],[2,4]]).transpose()
    bf = BlockFunction(bnd, f2, df2_dx)
    assert abs(bf.slab_vol(0.5,3.5)-2.75)<1e-12, 'Volume should be 2.75'
  
    # 1d function
    f3 = lambda x: x**2
    bnd = np.array([0,1])
    df3_dx = lambda x: 2*x
    bf = BlockFunction(bnd, f3, df3_dx)
    assert abs(bf.slab_vol(0,1)-0.75) < 1e-12
    assert abs(bf.slab_vol(0.5,1)-0.25) < 1e-12

def test_03():
    """
    One dimensional functions
    """
    f = lambda x: 2*x
    dfdx = lambda x: 2*np.ones(x.shape)
    n_intervals = 5
    bnd = np.array([0,2])
    g = GridFunction(bnd, n_intervals, f, dfdx)
    x_min, f_min, b_min, x_max, f_max, b_max = g.global_extrema()
    n_partition = 10
    f_grid = np.linspace(f_min,f_max,n_partition+1)
    f_prob = np.ones(n_partition)/n_partition
    g.set_output_histogram(f_grid, f_prob)
    g.compute_histogram()
    p = np.empty(n_intervals)
    count = 0
    for block in g.grid().values():
        p[count] = block.probability()
        count += 1
    assert abs(p.sum()-1)<1e-12, 'Probabilities should add up to 1.'    
    
    g.plot('function')
    g.plot('histogram')

if __name__ == '__main__':
    """
    
    """
    plt.close('all')
    test_01()
    #test_02()
    #test_03()
    """
    #
    # Define computational grid
    #
    x_min, x_max, y_min, y_max = 0, 2, 0, 2
    nx, ny = 20, 20  # number of cells in each direction
    x_grid = np.linspace(x_min,x_max,nx+1)
    y_grid = np.linspace(y_min,y_max,ny+1)
    bnd = np.array([[0,2],[0,2]])
    n_intervals = (nx, ny)
    #GridFunction(bnd, n_intervals, f, df_dx)
    #fig, ax = plt.subplots()
    #ax.xaxis.set_ticks(x_grid)
    #ax.yaxis.set_ticks(y_grid)    
    
    #ax.set_xlim([x_min, x_max])
    #ax.set_ylim([y_min, y_max])    
    #ax.grid(color='black')
    #plt.xticks(rotation=45)
    dx = (x_max-x_min)/nx
    dy = (y_max-y_min)/ny
    # Construct piecewise linear approximation of surface
    T = {}
    for j in range(nx):
        for i in range(ny):
            x_bnd = np.array([[x_min+i*dx,x_min+(i+1)*dx],
                              [y_min+j*dy,y_min+(j+1)*dy]])
            bij = BlockFunction(x_bnd, f, df_dx)
            bij.assign_address((i,j))
            T[(i,j)] = bij
    
    x_min, f_min, b_min, x_max, f_max, b_max = global_extrema(T)
    print(x_min, f_min)
    x_path_up, f_path_up = monotone_path(T, x_min, b_min, direction='up')
    x_path_down, f_path_down = monotone_path(T, x_max, b_max, direction='down')
    #t = transverse_parametrization(T) 
    #t = np.array(t)
    #plt.plot(t[:,0],t[:,1],'.-r')
    
    X, Y = np.meshgrid(x_grid, y_grid)
    levels = [0.06, 0.12, 0.18, 0.24, 0.3, 0.36, 0.42]
    xy = np.array([X.ravel(),Y.ravel()])
    cp = plt.contour(X,Y,f(xy).reshape(X.shape), levels)
    plt.clabel(cp, inline=True)
    for x in x_path_down:
        plt.plot(x[0,:], x[1,:], 'o-b')
        
    for x in x_path_up:
        plt.plot(x[0,:], x[1,:], '+-b')
    
    
    fig, ax = plt.subplots()
    count = 0
    for f in f_path_up:
        plt.plot(np.array([count,count+1]), f.squeeze())
        count += 1
        
    for f in f_path_down:
        plt.plot(np.array([count,count+1]), f.squeeze())
        count += 1
    #x_min, f_min, x_max, f_max = global_extrema(T)
    
    
    
    
    #
    # Estimate probability density function
    #
    # Assume a uniform distribution
    n_partitions = 10
    f_grid = np.linspace(f_min, f_max, n_partitions+1)
    PAj = np.ones(n_partitions)/n_partitions
    n_blocks = len(T)
    Aj = np.empty(n_partitions)
    V = np.empty((n_partitions,nx*ny))
    for j in range(n_partitions):
        b_Aj = np.empty((nx,ny))
        for i1 in range(nx):
            for i2 in range(ny):
                block = T[(i1,i2)]
                b_Aj[i1,i2] = block.slab_vol(f_grid[j],f_grid[j+1])
        Aj[j] = b_Aj.sum()
        V[j,:] = b_Aj.ravel()/Aj[j]    
    Pbi = np.dot(V.transpose(), PAj).reshape((nx,ny))
    print(Pbi.sum())
    fig  = plt.figure()
    ax = fig.gca(projection='3d')
    for i1 in range(nx):
        for i2 in range(ny):
            block = T[(i1,i2)]
            x = block.vertices()
            [X,Y] = np.meshgrid(x[0,:],x[1,:])
            Z = Pbi[i1,i2]*np.ones(X.shape)
            surf = ax.plot_surface(X,Y,Z, alpha=0.5, cmap='viridis', \
                                   vmin=Pbi.min(), vmax=Pbi.max())
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    for ix in range(nx):
        for iy in range(ny):
            block = T[(ix,iy)]
            x = block.vertices()
            [X,Y] = np.meshgrid(x[0,:],x[1,:])
            Z = block.evaluate(np.array([X.ravel(),Y.ravel()])).reshape(X.shape)
            surf = ax.plot_surface(X,Y,Z, cmap='viridis', vmin=f_min, vmax=f_max)
    
    plt.matshow(Pbi, cmap='Blues_r')    
       
       
    f = lambda x: np.exp(-x[0])
    dfdx = lambda x: -np.exp(-x[0])
    nx = 10
    x_min, x_max = 0, 3
    dx = (x_max-x_min)/nx
    Blocks = []
    for ix in range(nx):
        bnd = np.array([[x_min+ix*dx, x_min+(ix+1)*dx]])
        bi = BlockFunction(bnd, f, dfdx)
        bi.assign_address(ix)
        Blocks.append(bi)
    fig = plt.figure()
    for ix in range(nx):
        block = Blocks[ix]
        x_bnd = block.bnd()
        plt.plot(x_bnd, block.evaluate(x_bnd),'b')
    
    xm,ym = x_min + (ix+0.5)*dx, y_min + (iy+0.5)*dy
    q0 = qfn(xm,ym)
    dq_dx, dq_dy = grad_q(xm,ym)
    # Loop over vertices, evaluate linearization at each vertex             
    box_min = [x_min + ix*dx, y_min + iy*dy]
    box_max = [x_min + (ix+1)*dx, y_min + (iy+1)*dy]
    vertices = np.array([p for p in product(box_min, box_max)])
    for v in vertices:
        q = q0 + dq_dx*(v[0]-xm) + dq_dy*(v[1]-ym)
        print(q0)
        if q_max == None:
            q_max = q0
            xy_max = v
        elif q0 > q_max:
            q_max = q0
            xy_max = v
    
    plt.plot(xm,ym,'.b')            
    #
    # Compute local linearization at cell midpoint
    #             
    
    def f(x):
    '''
    Forward map from parameter space [0,2]x[0,2] to [0,oo].
    
    Inputs: 
    
        x: (2,k) vector of inputs
        
    Outputs: 
    
        fx: (k,) vector of outputs
    
    '''
    #
    # Check if dimensions match
    # 
    assert x.shape[0] == 2, \
        'Dimension of point differs from that of box.'
        
    # 
    # Turn into a column vector if necessary
    # 
    if len(x.shape) == 1:
        x = np.reshape(x, (-1,1))    
        
    #
    # Evaluate function
    #
    return x[0,:]*x[1,:]*np.exp(-(x[0,:]**2+1.25*x[1,:]**2-1))



def df_dx(x):
    '''
    Return the derivatives of f with respect to x, and y.
    
    Input:
    
        x: double, (2,k) array of k input vectors
        
    Outputs:
    
        dfdx: double (k,2) array of gradient vectors, one for each point
        
    '''
    #
    # Check if dimensions match
    # 
    assert x.shape[0] == 2, \
        'Dimension of point differs from that of box.'
        
    # 
    # Turn into a column vector if necessary
    # 
    if len(x.shape) == 1:
        x = np.reshape(x, (-1,1)) 
    
    #
    # Evaluate gradient
    # 
    x1, x2 = x[0,:], x[1,:]
    df_dx1 = x2*np.exp(-(x1**2+1.25*x2**2-1))*(1-2*x1**2)
    df_dx2 = x1*np.exp(-(x1**2+1.25*x2**2-1))*(1-2.5*x2**2)
    
    return np.array([df_dx1, df_dx2]).transpose().squeeze()
    
def transverse_path(T):
    '''
    Compute the transverse parametrization of a generalized contour
    
    Inputs: 
    
        T: BlockFunction, dictionary of hyperplanes that approximate the 
            function f within each subblock                        
            * keys: tuples specifying address of block in grid
            * values: BlockFunction object
            
    Output:
    
        t: double, (n,m) array of points defining the transverse curve
    '''
    path = []
    #
    # Compute the global minimum and -maximum
    # 
    x_min, f_min, block_min, x_max, f_max, block_max = global_extrema(T)
    
    #
    # Travel from global maximum along negative gradient
    #
    x_dwn, block_dwn = x_min, block_min
    path_dwn_x, path_dwn_f = monotone_path(T, x_dwn, block_dwn) 
    
    #
    # Travel from global max
    # 
    transverse_path = [x]
    count = 0
    while True:
        count += 1        
        
        p = -block.dfdx()  # descent direction
        x_bnd = block.bnd()
        
        #
        # Downwind boundary
        # 
        bnd_dwn = 0.5*(x_bnd[:,0]+x_bnd[:,1]) + \
                  0.5*np.sign(p)*(x_bnd[:,1]-x_bnd[:,0])
        
        
        #
        # Determine stepsize (how far to the nearest boundary?)
        #
        dx = np.squeeze(bnd_dwn-x)
        a_vec = np.where(p!=0, np.abs(dx/p), np.inf)
        i = np.argmin(a_vec)
        a = a_vec[i]
        
        #
        # 
        # 
        if np.abs(a) < 1e-13 or count > 20:
            break
        
        #
        # Update x
        # 
        x = x + a*p
        transverse_path.append(x)
        
        #
        # Update block
        #         
        address = np.array(block.address())
        address[i] += int(np.sign(p[i]))
        if any(address < 0):
            print('We have reached the boundary')
            break 
        print(tuple(address))
        block = T[tuple(address)]
        

    return transverse_path
        #new_address[i] += int(np.sign(p.squeeze()[i]))
        
 

def global_extrema(T):
    '''
    Compute the maximum/minimum of the piecewise linear approximation of 
    q over an n-dimensional rectangular grid.
    
    Input:
    
        T: BlockFunction, list of hyperplanes that approximate the function f
            within each sub-block. 
            * keys: tuples specifying address of block in grid
            * values: BlockFunction object
    Output:
    
        x_min, f_min, x_max, f_max: double, locations and function values of 
            global extrema.
    
    '''
    x_min, x_max, f_min, f_max = None, None, None, None
    b_min, b_max = None, None
    for block in T.values():
        #
        # Iterate over blocks
        #
    
        # Compute local extrema
        x_min_loc, f_min_loc, x_max_loc, f_max_loc = \
            block.extrema(min_or_max='both')
        
        # Update global minimum if necessary
        if (f_min is None) or f_min_loc < f_min:
            x_min, f_min, b_min  = x_min_loc, f_min_loc, block
        
        if (f_max is None) or f_max_loc > f_max:
            x_max, f_max, b_max = x_max_loc, f_max_loc, block
        
    return x_min, f_min, b_min, x_max, f_max, b_max
 

def monotone_path(T, x, block, direction='up'):
    '''
    Compute the monotone path transverse to the function's contours along 
    the path of steepest (i) increase (gradient), or (ii) decrease (-gradient). 
        
    Inputs:
    
        T: BlockFunction, dictionary of hyperplanes that approximate the 
            function f within each subblock. 
            * keys: tuples specifying address of block in grid
            * values: BlockFunction object
            
        x: x-value at which path begins
        
        block: BlockFunction, block in which initial point lies.
        
        direction: str, specifying whether the path goes along the gradient
            ('up') or against it ('down')
 
    
    Outputs: 
    
        path_x: double, list of (n,2) arrays containing the initial and final 
            x-values within each block. 
        
        path_f: double, list of (2,) arrays containing the approximate function 
            values at the points in 'path_x'
        
    '''  
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
        
        block = T[tuple(address)]
        x = x_new
        fx = fx_new        
        
        
       
"""
    
    
        


    
