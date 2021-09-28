#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  5 15:15:44 2021

@author: hans-werner
"""

import numpy as np
import Tasmanian

def example_10():

    print("\n---------------------------------------------------------------------------------------------------\n")
    print("Example 10: comparison between local polynomial and wavelet grids\n")

    iNumInputs = 2 # using two inputs for testing

    # test the error on a uniform dense grid with 10K points
    iTestGridSize = 100
    dx = np.linspace(-1.0, 1.0, iTestGridSize) # sample on a uniform grid
    aMeshX, aMeshY = np.meshgrid(dx, dx)
    aTestPoints = np.column_stack([aMeshX.reshape((iTestGridSize**2, 1)),
                                   aMeshY.reshape((iTestGridSize**2, 1))])

    def get_error(grid, model, aTestPoints):
        aGridResult = grid.evaluateBatch(aTestPoints)
        aModelResult = np.empty((aTestPoints.shape[0], 1), np.float64)
        for i in range(aTestPoints.shape[0]):
            aModelResult[i,:] = model(aTestPoints[i,:])
        return np.max(np.abs(aModelResult[:,0] - aGridResult[:,0]))

    def sharp_model(aX):
        return np.ones((1,)) * aX[0] / (1.0 + 100.0 * np.exp(-10.0 * aX[1]))

    grid_poly = Tasmanian.makeLocalPolynomialGrid(iNumInputs, 1, 5,
                                                  iOrder = 1, sRule = "localp")
    
    grid_wavelet = Tasmanian.makeWaveletGrid(iNumInputs, 1, 1, 1)

    print("            polynomial              wavelet")
    print("  points         error  points        error")

    fTolerance = 1.E-5;
    while((grid_poly.getNumNeeded() > 0) or (grid_wavelet.getNumNeeded() > 0)):
        Tasmanian.loadNeededPoints(lambda x, tid: sharp_model(x), grid_poly, 4);
        Tasmanian.loadNeededPoints(lambda x, tid: sharp_model(x), grid_wavelet, 4);

        # print the results at this stage
        print("{0:>8d}{1:>14.4e}{2:>8d}{3:>14.4e}".format(
            grid_poly.getNumLoaded(), get_error(grid_poly, sharp_model, aTestPoints),
            grid_wavelet.getNumLoaded(), get_error(grid_wavelet, sharp_model, aTestPoints)))

        # setting refinement for each grid
        grid_poly.setSurplusRefinement(fTolerance, -1, "fds");
        grid_wavelet.setSurplusRefinement(fTolerance, -1 , "fds");


if (__name__ == "__main__"):
    example_10()
