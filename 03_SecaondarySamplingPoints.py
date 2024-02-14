##############################################################################################
############### Probabilistic hazard map with secondary sampling points (SSPs) ###############
##############################################################################################
import numpy as np
import pandas as pd

from myFunctions import *

def main():
    ### command file name ###
    str_file = "pcq_template_v3.cmd"

    ### get array of parameters ###
    params_list = read_cmd(str_file)

    ### Read PC coefficients & make Probabilistic hazard map with uniform distributions
    load_coefficients(params_list, params_list[-1])

#   Read PC coefficients & make probabilistic hazard map with uniform distributions
##  input : p_list - array of parameters comes from read_cmd()
##          str_path - path to the directoly holding the file of quadrature points
# p_list =[N_val, x1_min, x1_max, ..., xD_min, xD_max, NP, NQ, NSSP, hcri, str_path]
def load_coefficients(p_list, str_path="./"):
    N_val = p_list[0]

    str_folder = str_path + "PCQ_bk/"
    NP, NQ = p_list[2*N_val+1], p_list[2*N_val+2]
    N_SSP, h_cri = p_list[2*N_val+3], p_list[2*N_val+4]
    str_path = p_list[2*N_val+5]

    trunc_K = NP + 1
    ### Set SSPs and evaluate them with Legendre polynomials (preparation)
    SSP_LIST = np.linspace(-1, 1, N_SSP)
    LEGENDRE_SSP_LIST = [special.eval_legendre(i*np.ones(N_SSP), SSP_LIST) for i in range(trunc_K)]

    ### Read csv file of b_k as matrix
    str_file = str_folder + "Bk_NP{:d}_NQ{:d}.csv".format(NP, NQ)
    coef = np.loadtxt(str_file, delimiter=",", skiprows=1)
    MAT_Bk = np.array(coef)

    ## Algebric calculation with SSP
    MAT_PHI = SSP_calc(N_val, trunc_K, LEGENDRE_SSP_LIST)

    ### Uncertainty output
    MAT_H = np.dot(MAT_Bk, MAT_PHI)

    ### Calculate exceed probability at each grid
    P_LIST = np.array([np.count_nonzero(h_array>h_cri) / N_SSP**N_val for h_array in MAT_H])

    ###save file
    header_dict = asc_header(str_path+"output_0.asc")
    str_file = str_path + "prob_NP{:d}-NQ{:d}.asc".format(NP, NQ)
    SaveAsc(np.transpose(P_LIST), str_file, header_dict)

#   Make tensor of SSPs
##  input : N_val  - Number of uncertainty inputs (less than 4)
##          trunc_K - Maximum degree of orthogonal polynomial
##          LEGENDRE_SSP_LIST - List of L_n(xi); xi in [-1, 1]
##              [[L_0(xi_1), L_0(xi_2), ..., L_0(xi_N_SSP)],, ..., [L_K(xi_1), L_K(xi_2), ..., L_K(xi_N_SSP)]]
def SSP_calc(N_val, trunc_K, LEGENDRE_SSP_LIST):
    if N_val==1:
        return np.array([
            [x1 for x1 in LEGENDRE_SSP_LIST[i1]]
            for i1 in range(trunc_K)
        ])
    elif N_val==2:
        return np.array([
            [x1*x2 for x1 in LEGENDRE_SSP_LIST[i1] for x2 in LEGENDRE_SSP_LIST[i2]]
            for i1 in range(trunc_K) for i2 in range(trunc_K-(i1))
        ])
    elif N_val==3:
        return np.array([
            [x1*x2*x3 for x1 in LEGENDRE_SSP_LIST[i1] for x2 in LEGENDRE_SSP_LIST[i2] for x3 in LEGENDRE_SSP_LIST[i3]]
            for i1 in range(trunc_K) for i2 in range(trunc_K-(i1)) for i3 in range(trunc_K-(i1+i2))
        ])
    else:
        print("Number of uncertainty inputs is assumed to be less than 3."); exit(1)

def SaveAsc(data1D, rec_file, header_dict):
    ### Set header
    x_num = header_dict["ncols"]
    y_num = header_dict["nrows"]
    asc_header = ""
    for mykey, myvalue in header_dict.items():
        asc_header += "{:s}\t\t{:.0f}\n".format(mykey, myvalue)

    ### 1D array to ascii format
    GRID_DATA = []
    for jj in range(y_num):
        str_rec = []
        for ii in range(x_num):
            i = jj*x_num + ii
            str_rec += ["{:.5f}".format(data1D[i])]
        GRID_DATA += [str_rec]
    ### Save probability as asc file
    with open(rec_file, "w") as f:
        f.write(asc_header)
        for l in GRID_DATA:
            for ll in l:
                f.write(ll); f.write(" ")
            f.write("\n")

if __name__ == "__main__":
    main()

exit(1)
