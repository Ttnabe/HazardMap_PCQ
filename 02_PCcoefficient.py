"""
Calculate the PC coefficients for 2 variables and hazard map with uniform SSPs
Created on Oct. 12th, 2023, by Takahiro Tanabe
"""
import numpy as np
import pandas as pd
from myFunctions import *

def main():
    str_file = "pcq_template_v2.cmd"

    ### Get parameter list
    params_list = read_cmd(str_file)

    ### Calculate the coefficients of PCE
    CalcCoef(params_list)

##  Calculate PC coefficients & save them to csv file
##  input : p_list - array of parameters comes from read_cmd()
# p_list =[N_val, x1_min, x1_max, ..., xD_min, xD_max, NP, NQ, NSSP, cri, str_path, rec_path]
def CalcCoef(p_list):
    ### Number of uncertainty inputs
    N_val = p_list[0]

    rec_path, str_path = p_list[-1], p_list[-2]
    rec_folder = rec_path + "PCQ_bk/"
    os.makedirs(rec_folder, exist_ok=True)

    x_sum, x_diff = [], []
    for n in range(N_val):
        x_sum += [p_list[2*(n+1)] + p_list[2*n+1]]
        x_diff += [p_list[2*(n+1)] - p_list[2*n+1]]
    NP, NQ = p_list[2*N_val+1], p_list[2*N_val+2]
    #MODIFY--------------------------
    #header = asc_header(str_path+"output_0.asc")
    header = asc_header(str_path+"../10_RES/0-50_h.asc")
    #--------------------------MODIFY
    N_DATA = header["ncols"] * header["nrows"]

    ## csv file recording nodes and weight (made by 01_Set*.py)
    str_file = rec_path + "{:d}param_PCQ{:d}.csv".format(N_val, NQ)
    INPUT_LIST = np.loadtxt(str_file, delimiter=",", skiprows=1)

    # Calculate PC coefficients
    ### Array whose component is [[xi_1^1, ..., xi_D^1], ..., [xi_1^NQ, ..., xi_D^NQ]] for all the quadrature points (x1, ..., xD) (size: (NQ^N_val, N_val))
    INPUT_LIST_XI = np.transpose((2*np.transpose(INPUT_LIST)[0::2] - np.array([np.full(len(INPUT_LIST), x) for x in x_sum])) / np.array([np.full(len(INPUT_LIST), x) for x in x_diff]))
    ### Array whose component is [[w_1^1,...,w_D^1], ..., [w_1^NQ, w_D^NQ]] for all weight (w1, ..., wNQ) (size: (NQ^N_val, N_val))
    WEIGHT_LIST = np.array([x[1::2] for x in INPUT_LIST])

    trunc_K = NP + 1##"trunc_K" is the upper limit of polynomials
    ### Array whose component is
    ### [[(1/LAMBDA)*(w_1^1*Le_0(xi_1^1))*...*(w_D^1*Le_0(xi_D^1)), ..., (1/LAMBDA)*(w_1^NQ*Le_0(xi_1^NQ))*...*(w_D^NQ*Le_0(xi_D^NQ))],
    ###  ...
    ###  [(1/LAMBDA)*(w_1^1*Le_K(xi_1^1))*...*(w_D^1*Le_K(xi_D^1)), ..., (1/LAMBDA)*(w_1^NQ*Le_K(xi_1^NQ))*...*(w_D^NQ*Le_K(xi_D^NQ))]]
    ### with K: trunc_K, Le_n: n-th Legendre polynomial (size: (K, NQ^N_val))
    if N_val==1:
        LEGENDRE_XI_LIST = np.array([
                            [(1./LAMBDA_n(i1)) * (q*legendre(i1, x)) for (x, q) in zip(INPUT_LIST_XI, WEIGHT_LIST)]
                            for i1 in range(trunc_K)])
    elif N_val==2:
        LEGENDRE_XI_LIST = np.array([
                            [(1./(LAMBDA_n(i1)*LAMBDA_n(i2))) * (q[0]*legendre(i1, x[0])) * (q[1]*legendre(i2, x[1])) for (x, q) in zip(INPUT_LIST_XI, WEIGHT_LIST)]
                            for i1 in range(trunc_K) for i2 in range(trunc_K-i1)])
    elif N_val==3:
        LEGENDRE_XI_LIST = np.array([
                            [(1./(LAMBDA_n(i1)*LAMBDA_n(i2)*LAMBDA_n(i3))) * (q[0]*legendre(i1, x[0])) * (q[1]*legendre(i2, x[1])) * (q[2]*legendre(i3, x[2])) for (x, q) in zip(INPUT_LIST_XI, WEIGHT_LIST)]
                            for i1 in range(trunc_K) for i2 in range(trunc_K-i1) for i3 in range(trunc_K-(i1+i2))])
    else:
        print("Number of uncertainty inputs is assumed to be less than 3."); exit(1)

    ### Column name "b_k"
    str_col = ["b_{:d}".format(i) for i in range(len(LEGENDRE_XI_LIST))]

    ### Array of output (H_max in Appendix A of Tanabe et al.):
    ### [[output(x0), output(x1), ..., output(xM)]_\theta_1,]
    ###  ....
    ###  [output(x0), output(x1), ..., output(xM)]_\theta_(NQ^N_val)]
    ### x0, x1, ..., xM indicate cell position
    ### \theta_d (d=1, ..., NQ^N_val) indicates input variable
    #MODIFY--------------------------
    #MATRIX_OUT = np.transpose([np.ravel(np.loadtxt(str_path+"output_{:d}.asc".format(nq), skiprows=6)) for nq in range(len(INPUT_LIST))])
    MATRIX_OUT = np.transpose([np.ravel(np.loadtxt(str_path+"csv_Hmax{:d}.csv".format(nq), delimiter=",")) for nq in range(len(INPUT_LIST))])
    #--------------------------MODIFY
    ### Calculate PC coefficients (Eq.(A2) in Appendix A of Tanabe et al.):
    ### [[b_0(x0), b_1(x0), ..., b_K(x0)],...,[b_0(xM), b_1(xM), ..., b_K(xM)]]
    BK_LIST = [[np.inner(out, phi_k) for phi_k in LEGENDRE_XI_LIST] for out in MATRIX_OUT]

    ## Save PC coefficients as csv file
    ## col: the order of PC coefficient, row: position
    df = pd.DataFrame(BK_LIST, columns=str_col)
    rec_file = rec_folder + "Bk_NP{:d}_NQ{:d}.csv".format(NP, NQ)
    df.to_csv(rec_file, index=False)

##  Calculate int{Le_n*Le_n}
##  input : n - degree of Legendre polynomial
def LAMBDA_n(n):
    return 2. / (2*n + 1)


if __name__ == "__main__":
    main()

exit(1)
