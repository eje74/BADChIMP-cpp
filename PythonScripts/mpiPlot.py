import numpy as np
import matplotlib.pyplot as plt

num_iter = 100
num_proc = 3
sys_size = (7, 12)

rho0 = np.zeros(sys_size);
rho1 = np.zeros(sys_size);

input_dir = "/home/ejette/Programs/GITHUB/badchimpp/input/"

for iter in np.arange(num_iter):
    for proc in np.arange(num_proc):
        file_name = "rho_val_" + str(proc) + "_" + str(iter) + ".dat"
        dta = np.loadtxt(input_dir + file_name)
        for n in np.arange(dta.shape[0]):
            rho0[int(dta[n,1]), int(dta[n,0])] = dta[n, 2]
            rho1[int(dta[n,1]), int(dta[n,0])] = dta[n, 3]

        fig = plt.figure(0)
        plt.pcolormesh(rho0)
        fig.savefig(input_dir + "rho0_"+str(iter)+".png", format='png')

        fig = plt.figure(1)
        plt.pcolormesh(rho1)
        fig.savefig(input_dir + "rho1_"+str(iter)+".png", format='png')
