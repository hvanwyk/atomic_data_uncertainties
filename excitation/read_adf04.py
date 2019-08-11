#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 22:34:39 2019

@author: kyle
"""

import numpy as np
import pandas as pd
import re

def read_adf04(filename):
    with open(filename, "r") as adf04:
        adf04.readline()
        levels = []
        line = adf04.readline()
        while "-1" not in line:
            line = line.strip()
            ang = re.findall("\(\d\)\d\( \d.\d\)", line)
            line = re.sub("\(\d\)\d\( \d.\d\)", "", line)
            conf = re.findall("(\d[a-zA-Z]\d)", line)
            conf = " ".join(conf)
            line = re.sub("(\d[a-zA-Z]\d)", "", line)
            line = line.split()
            
            
            current = [line[0]]
            current.append(conf)
            current += ang
            current.append(line[1])
            levels.append(current)
            
            line = adf04.readline()
        T = adf04.readline().split()
        line = adf04.readline()
        r = re.findall("\-?\d.\d{2}[\+|\-]\d{2}", line)
        line = re.sub("\-?\d.\d{2}[\+|\-]\d{2}", "", line)
        arr = line.split()
        arr += r
        
        rates = np.array(arr).reshape((1,-1))

        line = adf04.readline()
        while line.strip() != "-1":
            r = re.findall("\-?\d.\d{2}[\+|\-]\d{2}", line)
            line = re.sub("\-?\d.\d{2}[\+|\-]\d{2}", "", line)
            arr = line.split()
            arr += r
            rates = np.concatenate((rates, np.array(arr).reshape((1,-1))),axis=0)
            line = adf04.readline()
            
    levels = np.array(levels)
    levels = pd.DataFrame({"#": levels[:, 0], "config": levels[:, 1], "(2S+1)L( 2J)": levels[:, 2], "Energy": levels[:, 3]})
    T = np.array(T)
    dct = {"initial": rates[:, 0], "final": rates[:, 1]}
    for i in range(2, rates.shape[1]):
        dct[f"T{i-2}"] = rates[:, i]

    rates = pd.DataFrame(dct)
    
    return levels, T, rates

filename = "isoelectronic/he-like/o6/adas/adf04"

levels, T, rates = read_adf04(filename)

print(levels)


