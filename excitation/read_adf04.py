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
        hdr = adf04.readline()
        levels = []
        line = adf04.readline()
        while "-1" not in line:
            hdr += line
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
        hdr += line
        
        line = adf04.readline()
        hdr += line
        T = line.split()
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
    
    shape = rates.shape

    rates = rates.ravel()
    for i in range(len(rates)):
        if "+" in rates[i]:
            r = rates[i].split("+")
            r = float("e+".join(r))
            rates[i] = r
        elif "-" in rates[i]:
            neg = False
            if rates[i][0] == "-":
                neg = True
            r = rates[i].split("-")
            if neg:
                r = float("e-".join(r[1:]))
                r *= -1
            else:
                r = float("e-".join(r))
            rates[i] = r
    for i in range(len(T)):
        if "+" in T[i]:
            t = T[i].split("+")
            t = float("e+".join(t))
            T[i] = t
        elif "-" in T[i]:
            neg = False
            if T[i][0] == "-":
                neg = True
            t = T[i].split("-")
            if neg:
                t = "e-".join(t[1:])
                t *= -1
            else:
                t = "e-".join(t)
            T[i] = t
        else:
            T[i] = float(T[i])
    
    rates = rates.reshape(shape)
    levels = np.array(levels)
    T = np.array(T)
    levels = pd.DataFrame({"#": levels[:, 0], "config": levels[:, 1], "(2S+1)L( 2J)": levels[:, 2], "Energy": levels[:, 3]})

    return levels, T.astype(np.float)[2:], rates.astype(np.float), hdr

def rates_dataframe(filename):
    
    levels, T, rates, hdr = read_adf04(filename)
    
    dct = {"initial": rates[:, 0], "final": rates[:, 1]}
    dct["A"] = rates[:, 2]
    for i in range(3, rates.shape[1]-1):
        dct[f"{T[i-3]}"] = rates[:, i]
    dct["inf_energy"] = rates[:,-1]

    rates_df = pd.DataFrame(dct)
    rates_df.sort_values(by=["final", "initial"])
    
    return rates_df
    
def compare(file_1, file_2):
    lev1, T1, r1 = read_adf04(file_1)
    lev2, T2, r2 = read_adf04(file_2)
    
    r1 = r1.ravel()
    r2 = r2.ravel()
    
    r1[r1==0] = 1e-99
    r2[r2==0] = 1e-99
    
    diff = np.abs((r1-r2)/((r1+r2)/2))
    
    
    max_diff = np.max(diff)
    print(f"% difference: {max_diff*100}")

    return diff

def compare_ground(file_1, file_2):
    df_1 = rates_dataframe(file_1)
    df_1 = df_1.loc[df_1["final"]==1.0]
    df_2 = rates_dataframe(file_2)
    df_2 = df_2.loc[df_2["final"]==1.0]

    r_1 = df_1.values[:,8:]
    r_2 = df_2.values[:,8:]
    r_1[r_1==0] = 1e-99
    r_2[r_2==0] = 1e-99
    
    diff = np.abs((r_1-r_2)/((r_1 + r_2)/2))
    
    max_diff = np.max(diff)
    
    print(f"% difference for up to n=4 transitions to ground: {max_diff*100}")
    return diff
    
    
    
    
if __name__ == "__main__":
    
    file_1 = "adf04"
    
    data = read_adf04(file_1)
    
    df = rates_dataframe(file_1)
    df = df[df["final"]==1]
    print(df)
    print(data[1][9])
    
    
