import numpy as np
import matplotlib.pyplot as plt

num_iter = 20001
write_interval = 100
start_iter = 0

num_proc = 3
sys_size = (8, 200, 100)

rho0 = np.zeros(sys_size);
rho1 = np.zeros(sys_size);

input_dir = "/home/ejette/Programs/GITHUB/badchimpp/output/"

for iter in np.arange(start_iter, num_iter,  write_interval):

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
    if iter == 0:
        rhoSolid = np.copy(rho0[4, :, :])

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

    fig = plt.figure(0)
    plt.clf()
#    plt.pcolormesh(rho0[4, : , :])
    c0 = rho0[4, :, :] / (rho0[4, :, :] + rho1[4, :, :])
    maskSolid = np.ma.masked_where(rhoSolid < 0.1,c0)
    plt.pcolormesh(maskSolid)
    plt.colorbar()
    plt.axis('tight')
    plt.axis('scaled')
    fig.savefig(input_dir + "c0_"+str(iter)+".png", format='png')
