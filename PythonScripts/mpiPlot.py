import numpy as np
import matplotlib.pyplot as plt

begin_iter = 67600
num_iter = 200001
write_interval = 200
num_proc = 3
sys_size = (8, 200, 100)

rho0 = np.zeros(sys_size);
rho1 = np.zeros(sys_size);
q0 = np.zeros(sys_size); 
q1 = np.zeros(sys_size);

input_dir = "/home/ejette/Programs/GitHub/BADChIMP-cpp/output/"

for proc in np.arange(num_proc):
    file_name = "rho_val_" + str(proc) + "_" + str(0) + ".dat"
rhoBackGround = np.ones(rho0.shape)
rhoBackGround[rho0 < 0.1] = 0

for iter in np.arange(begin_iter, num_iter,  write_interval):
    for proc in np.arange(num_proc):
        file_name = "rho_val_" + str(proc) + "_" + str(iter) + ".dat"
        dta = np.loadtxt(input_dir + file_name)
        # for n in np.arange(dta.shape[0]):
        #     rho0[int(dta[n,2]), int(dta[n,1]), int(dta[n,0])] = dta[n, 3]
        #     rho1[int(dta[n,2]), int(dta[n,1]), int(dta[n,0])] = dta[n, 4]
        x = dta[:,0].astype(int)
        y = dta[:,1].astype(int)
        z = dta[:,2].astype(int)
        rho0[z, y, x] = dta[:, 3]
        rho1[z, y, x] = dta[:, 4]
        q0[z, y, x] = dta[:, 5]
        q1[z, y, x] = dta[:, 6]

    fig = plt.figure(0)
    plt.clf()
#    plt.pcolormesh(rho0[4, : , :])
    plt.pcolormesh(np.sum(rho0, axis=0))
    plt.colorbar()
    plt.axis('tight')
    plt.axis('scaled')
    fig.savefig(input_dir + "rho0_"+str(iter)+".png", format='png')

    fig = plt.figure(1)
    plt.clf()
#    plt.pcolormesh(rho1[4, : , :])
    plt.pcolormesh(np.sum(rho1, axis=0))
    plt.colorbar()
    plt.axis('tight')
    plt.axis('scaled')
    fig.savefig(input_dir + "rho1_"+str(iter)+".png", format='png')

    fig = plt.figure(2)
    plt.clf()
#    plt.pcolormesh(rho1[4, : , :])
    plt.pcolormesh(np.sum(q0, axis=0))
    plt.colorbar()
    plt.axis('tight')
    plt.axis('scaled')
    fig.savefig(input_dir + "q0_"+str(iter)+".png", format='png')

    fig = plt.figure(3)
    plt.clf()
#    plt.pcolormesh(rho1[4, : , :])
    plt.pcolormesh(np.sum(q1, axis=0))
    plt.colorbar()
    plt.axis('tight')
    plt.axis('scaled')
    fig.savefig(input_dir + "q1_"+str(iter)+".png", format='png')

    fig = plt.figure(4)
    plt.clf()
#    plt.pcolormesh(rho0[4, : , :])
    plt.pcolormesh(np.sum(rho0 + rho1 - rhoBackGround, axis=0))
    plt.colorbar()
    plt.axis('tight')
    plt.axis('scaled')
    fig.savefig(input_dir + "rhoDiff_"+str(iter)+".png", format='png')
