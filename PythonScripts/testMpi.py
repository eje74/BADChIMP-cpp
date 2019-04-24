import numpy as np
import matplotlib.pyplot as plt

num_iter = 1001
num_proc = 1
sys_size = (8, 800, 400)

rho0_n3 = np.zeros(sys_size);
rho1_n3 = np.zeros(sys_size);

rho0_n1 = np.zeros(sys_size);
rho1_n1 = np.zeros(sys_size);


input_dir_n3 = "/home/ejette/Programs/GITHUB/badchimpp/output/N3/"
#input_dir_n1 = "/home/ejette/Programs/GITHUB/badchimpp/output/N1/"
input_dir_n1 = "/home/ejette/Programs/GITHUB/badchimpp/output/"


for iter in [num_iter -1]:#np.arange(num_iter):
    # # nproc = 3
    # for proc in np.arange(num_proc):
    #     file_name = "rho_val_" + str(proc) + "_" + str(iter) + ".dat"
    #     dta = np.loadtxt(input_dir_n3 + file_name)
    #     for n in np.arange(dta.shape[0]):
    #         rho0_n3[int(dta[n,1]), int(dta[n,0])] = dta[n, 2]
    #         rho1_n3[int(dta[n,1]), int(dta[n,0])] = dta[n, 3]

    # nproc = 1
    file_name = "rho_val_" + str(0) + "_" + str(iter) + ".dat"
    dta = np.loadtxt(input_dir_n1 + file_name)
    for n in np.arange(dta.shape[0]):
        rho0_n1[int(dta[n,2]), int(dta[n,1]), int(dta[n,0])] = dta[n, 3]
        rho1_n1[int(dta[n,2]), int(dta[n,1]), int(dta[n,0])] = dta[n, 4]

#    print(max(np.abs(rho0_n3.flatten() - rho0_n1.flatten()) ) )

#np.set_printoptions(precision=3)
#print(rho0_n1)
#print(rho1_n1)
