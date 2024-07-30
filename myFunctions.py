"""
Module: Common functions among PCQ procedures
Created on Oct. 12th, 2023, by Takahiro Tanabe
"""
import os
import sys
from scipy import special
import numpy as np

# Read command file & return the parameters list
## input : str_file - command file name
## output: list of parameters written in command file
def read_cmd(str_file):
    cmd_list = []
    ## Set parameter ranges & number of nodes for the quadature
    ## Read command file
    with open(str_file) as f:
        l = f.readlines()
        ## Number of uncertain variables
        N_val = int(l[0].split()[-1])
        cmd_list += [N_val]
        ### Minimum & maximum bounds of each parameter
        for n in range(N_val):
            cmd_list += [float(l[3+2*n].split()[-1])]
            cmd_list += [float(l[3+(2*n+1)].split()[-1])]

        ### N_P, N_Q, N_SSP
        for n in range(6+2*N_val, (6+2*N_val)+3):
            cmd_list += [int(l[n].split()[-1])]

        ### Threshold, folder path for output asc file (numerical result) & folder path for save
        cmd_list += [float(l[11+2*N_val].split()[-1])]
        cmd_list += [l[14+2*N_val][:-1]]
        cmd_list += [l[15+2*N_val][:-1]]
    return cmd_list


#   Get dictionary of asc-header
##  input : str_file - path for asc file
##  output: Dictionary-type variable whose components is header of asc file
def asc_header(str_file):
    header_dict = {}
    with open(str_file) as f:
        lines = f.read()
        for l in lines.split("\n")[:6]:
            blocks = l.split(" ")
            if blocks[0][0]=="n":
                header_dict[blocks[0]] = int(blocks[-1])
            else:
                header_dict[blocks[0]] = round(float(blocks[-1]), 2)
    return header_dict


##  Calculate Le_n(x)
##  input : n - degree of Legendre polynomial
##          x - A variable in [-1, 1]
def legendre(n, x):
    return special.eval_legendre(n, x)
