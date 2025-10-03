import numpy as np
import h5py
import matplotlib.pyplot as plt
import os

dir = '/Users/formanet/github/MCHPC/'

def evaluate_ampere(filename):
    fields = h5py.File(filename) 
    
    bz = fields["bz"][:]
    jy = fields["jy"][:]
    x  = fields["x"][:]
    
    fields.close()
    
    # --- analytic prediction ---
    x_start = x[0]
    x_end   = x[-1]
    L = x_end - x_start
    jy_pred = -(2*np.pi / L) * np.cos(2*np.pi*(x - x_start)/L)
    
    # --- plotting ---
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10,8), sharex=True)

    ax1.plot(x, bz, color='blue', linewidth=2)
    ax1.set_ylabel(r'$B_z$')
    ax1.set_title('Result of Ampere test')
    ax1.grid(True)
    
    ax2.plot(x, jy, color='red', label='Numerical $j_y$')
    ax2.plot(x, jy_pred, color='black', linestyle='--', label='Analytic $j_y$')
    ax2.set_xlabel(r'$x$')
    ax2.set_ylabel(r'$j_y$')
    ax2.legend()
    ax2.grid(True)
    
    plt.tight_layout()
    if not os.path.exists(os.path.join(dir,'IRT/plots/')):
        os.makedirs(os.path.join(dir,'IRT/plots/'))
    plt.savefig(os.path.join(dir,'IRT/plots/ampere_test.png'))
    plt.show()


def run_test():
    test_subdir = 'build/tests/ampere_faraday'
    evaluate_ampere(os.path.join(dir, test_subdir, 'sin_bz.h5'))
    

if __name__ == "__main__":
    run_test()