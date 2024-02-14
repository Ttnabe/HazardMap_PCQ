###############################################
#### Common functions among PCQ procedures ####
###############################################
import os
import sys
import numpy as np

# read command file & return the parameters written in the file
## input : str_file - name of command file
##         N_val - number of uncertainty input variables
## output: array of parameters written in command file
def read_cmd(str_file):
    cmd_list = []
    ## set parameter range & number of nodes for the quadature
    with open(str_file) as f:
        l = f.readlines()
        N_val = int(l[0].split()[-1])
        cmd_list += [N_val]
        ### minimum & maximum for limit of parameter interval
        for n in range(N_val):
            #print(n, l[3+2*n].split()[-1])
            cmd_list += [float(l[3+2*n].split()[-1])]
            #print(n, l[3+2*n+1].split()[-1])
            cmd_list += [float(l[3+(2*n+1)].split()[-1])]

        ### N_P, N_Q, N_SSP
        for n in range(6+2*N_val, 6+(2*N_val+3)):
            cmd_list += [int(l[n].split()[-1])]

        ### h_cri, path for output asc file (numerical result)
        cmd_list += [float(l[11+2*N_val].split()[-1])]
        cmd_list += [l[14+2*N_val][:-1]]
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
