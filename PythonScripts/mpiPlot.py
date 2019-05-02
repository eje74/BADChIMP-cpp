import numpy as np
import matplotlib.pyplot as plt

num_iter = 500001
write_interval = 1000
start_iter = 62000
num_proc = 3
sys_size = (8, 800, 400)

rho0 = np.zeros(sys_size);
rho1 = np.zeros(sys_size);

input_dir = "/home/olau/Programs/Git/BADChIMP-cpp/output/"

for iter in np.arange(start_iter, num_iter,  write_interval):
    
    for proc in np.arange(num_proc):
        file_name = "rho_val_" + str(proc) + "_" + str(iter) + ".dat"
        dta = np.loadtxt(input_dir + file_name)
        for n in np.arange(dta.shape[0]):
            rho0[int(dta[n,2]), int(dta[n,1]), int(dta[n,0])] = dta[n, 3]
            rho1[int(dta[n,2]), int(dta[n,1]), int(dta[n,0])] = dta[n, 4]

    if iter==0:
        rho0_0=rho0


    fig = plt.figure(0,figsize=(16,8))
    ax = fig.add_subplot(111)
    plt.clf()
    plt.pcolormesh(rho0[4, : , :])
    #plt.colorbar()
    plt.axis('tight')
    #ax.set_aspect('equal')
    plt.axis('scaled')
    fig.savefig(input_dir + "rho0_"+str(iter)+".png", format='png', bbox_inches='tight')

    fig = plt.figure(1,figsize=(16,8))
    ax = fig.add_subplot(111)
    plt.clf()
    plt.pcolormesh(rho1[4, : , :])
    #plt.colorbar()
    plt.axis('tight')
    #ax.set_aspect('equal')
    plt.axis('scaled')
    fig.savefig(input_dir + "rho1_"+str(iter)+".png", format='png', bbox_inches='tight')

    print(iter)