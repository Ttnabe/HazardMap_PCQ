#####################################################################
####### Set quadrature points (nodes) for Gaussian quadrature #######
#####################################################################
import itertools
from myFunctions import *

def main():
    ### command file name ###
    str_file = "pcq_template_v3.cmd"

    ### get array of parameters ###
    params_list = read_cmd(str_file)

    ### save quadrature points & weights ###
    save_nodes(params_list)

##  save quadrature points & weights
##  input : p_list - array of parameters comes from read_cmd()
##          N_val  - Number of uncertainty inputs
# p_list =[x1_min, x1_max, x2_min, x2_max, x3_min, x3_max,
#         NP, NQ, NSSP, hcri, str_path]
def save_nodes(p_list):
    ### Number of uncertainty inputs
    N_val = p_list[0]

    ### Calculate sum. & diff. of parameter edges for preparation
    x_sum, x_diff = [], []
    for n in range(N_val):
        x_sum += [(p_list[2*(n+1)]+p_list[2*n+1])/2]
        x_diff += [(p_list[2*(n+1)]-p_list[2*n+1])/2]
    NQ = p_list[2*(N_val+1)]; str_path = p_list[2*N_val+5]
    os.makedirs(str_path, exist_ok=True)
    print(x_sum, x_diff)
    ### calculate zero-points & each weight
    zl, wl = leg_weights_roots(NQ)
    ### convert xi\in[-1, 1] into parameter space x\in[X_min, X_max]
    y = np.array(zl)*np.array([np.full(len(zl), x) for x in x_diff]) + np.array([np.full(len(zl), x) for x in x_sum])
    ### header for csv file
    str_header = "".join(["param {0:d},weight {0:d},".format(n) for n in range(N_val)])[:-1] + "\n"

    ### iteration list
    i_list = list(itertools.product(range(NQ), repeat=N_val))
    ### save input variables to csv
    str_file = str_path + "{:d}param_PCQ{:d}.csv".format(N_val, NQ)
    with open(str_file, "w") as f:
        f.write(str_header)

        for i in i_list:
            str_f = "".join(["{:.7f},{:.7f},".format(yy[ii], wl[ii]) for ii, yy in zip(i, y)])[:-1] + "\n"
            f.write(str_f)

if __name__ == "__main__":
    main()
