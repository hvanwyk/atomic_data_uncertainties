import unittest
from density_of_inverse import BlockFunction, GridFunction
import numpy as np

class TestDensityOfInverse(unittest.TestCase):
    """
    """
    def test01(self):
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
        
        
        
    def test_02(self):
        """
        Test the computation of volumes of slabs
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
    
    
    def test_03(self):
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
    
if __name__ == "__main__":
    unittest.main()