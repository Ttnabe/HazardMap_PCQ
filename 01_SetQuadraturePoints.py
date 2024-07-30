"""
Set quadrature points (nodes) for Gaussian quadrature
Created on Oct. 12th, 2023, by Takahiro Tanabe
"""
import itertools
import sympy as sy
from myFunctions import *

def main():
    ### Command file name with path ###
    str_file = "pcq_template_v2.cmd"

    ### Get array of parameters ###
    params_list = read_cmd(str_file)

    ### Save quadrature points & weights ###
    save_nodes(params_list)

##  Save quadrature points & weights
##  input : p_list - array of parameters comes from read_cmd()
# p_list =[N_val, x1_min, x1_max, x2_min, x2_max, x3_min, x3_max,
#         N_P, N_Q, N_SSP, cri, str_path, rec_path]
def save_nodes(p_list):
    ### Number of uncertainty inputs
    N_val = p_list[0]

    ### Calculate sum. & diff. of parameter edges for preparation
    x_sum, x_diff = [], []
    for n in range(N_val):
        x_sum += [0.5*(p_list[2*(n+1)]+p_list[2*n+1])]
        x_diff += [0.5*(p_list[2*(n+1)]-p_list[2*n+1])]
    NQ = p_list[2*(N_val+1)]

    rec_path = p_list[-1]
    os.makedirs(rec_path, exist_ok=True)
    print(x_sum, x_diff)
    ### Calculate zero-points & each weight
    zl, wl = leg_weights_roots(NQ)
    ### Convert xi\in[-1, 1] into parameter space x\in[X_min, X_max]
    y = np.array(zl)*np.array([np.full(len(zl), x) for x in x_diff]) + np.array([np.full(len(zl), x) for x in x_sum])
    ### Header for csv file
    str_header = "".join(["param {0:d},weight {0:d},".format(n) for n in range(N_val)])[:-1] + "\n"

    ### Iteration list
    i_list = list(itertools.product(range(NQ), repeat=N_val))
    ### Save input variables to csv
    rec_file = rec_path + "{:d}param_PCQ{:d}.csv".format(N_val, NQ)
    with open(rec_file, "w") as f:
        f.write(str_header)

        for i in i_list:
            str_f = "".join(["{:.7f},{:.7f},".format(yy[ii], wl[ii]) for ii, yy in zip(i, y)])[:-1] + "\n"
            f.write(str_f)

# Calculate zero-points x_i of n-th order Legendre polynomial & its weight w_i for Gaussian quadature
## input : n - number of points (corresponding to NQ)
## output: arrays of zero points x_i=[x_1,...,x_n] & weight w_i=[w_1,...,w_n]
def leg_weights_roots(n):
    x = sy.Symbol('x')
    ## roots of n-th order Legendre polynomial
    roots = sy.Poly(sy.legendre(n, x)).all_roots()

    ## get x_i & its weight as array
    x_i = [rt.evalf(20) for rt in roots]
    w_i = [(2*(1-rt**2)/(n*sy.legendre(n-1, rt))**2).evalf(20) for rt in roots]
    return x_i, w_i

if __name__ == "__main__":
    main()
