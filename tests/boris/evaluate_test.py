import numpy as np
import h5py
import matplotlib.pyplot as plt
import os

dir = '/Users/formanet/github/MCHPC/'

def evaluate_test(filename, xideal = 5.05, videal = (0., 2., 0.), no_plot = False):
    uniform_bz = h5py.File(filename) 

    x  = uniform_bz["x"][:]
    vx = uniform_bz["vx"][:]
    vy = uniform_bz["vy"][:]
    vz = uniform_bz["vz"][:]
    uniform_bz.close()
    nstep = np.arange(len(x))

    # compare the final value with the analytical result
    delta_x = x[-1] - xideal
    delta_vx = vx[-1] - videal[0]
    delta_vy = vy[-1] - videal[1]
    delta_vz = vz[-1] - videal[2]
    # print results and assess numerical accuracy
    threshold = 1e-2 # arbitrary choice
    print("-" * 70)
    print(f"File: {filename}")
    print(f"Final position x: {x[-1]:.6f}, Ideal: {xideal:.6f}, Delta: {delta_x:.6e}")
    print(f"Final velocity vx: {vx[-1]:.6f}, Ideal: {videal[0]:.6f}, Delta: {delta_vx:.6e}")
    print(f"Final velocity vy: {vy[-1]:.6f}, Ideal: {videal[1]:.6f}, Delta: {delta_vy:.6e}")
    print(f"Final velocity vz: {vz[-1]:.6f}, Ideal: {videal[2]:.6f}, Delta: {delta_vz:.6e}")
    if (abs(delta_x) < threshold and 
        abs(delta_vx) < threshold and 
        abs(delta_vy) < threshold and 
        abs(delta_vz) < threshold):
        print("Test PASSED: Numerical solution is within the acceptable threshold.")
    else:
        print("Test FAILED: Numerical solution exceeds the acceptable threshold.")
    print("-" * 70)

    if not no_plot:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10,8), sharex=True)

        # Position vs time
        ax1.plot(nstep, x, color='blue', linewidth=2)
        ax1.set_ylabel("Position x")
        ax1.set_title("1D Particle Motion with 3D Velocity")
        ax1.grid(True)

        # Velocity vs time
        ax2.plot(nstep, vx, label='vx', color='red')
        ax2.plot(nstep, vy, label='vy', color='green')
        ax2.plot(nstep, vz, label='vz', color='purple')
        ax2.set_xlabel("Time step")
        ax2.set_ylabel("Velocity")
        ax2.legend()
        ax2.grid(True)

        plt.tight_layout()
        # save Figure 
        if not os.path.exists(os.path.join(dir,'IRT/plots/')):
            os.makedirs(os.path.join(dir,'IRT/plots/'))
        
        figname = filename.split('/')[-1].replace('.h5', '.png')
        plt.savefig(os.path.join(dir,'IRT/plots/',figname))
        # plt.close()

def run_test():
    # Call the function for each file
    # uniform field - test particle completes two cycles and returns to its initial state 
    test_subdir = 'build/tests/boris'
    evaluate_test(os.path.join(dir, test_subdir, 'uniform_bz.h5'), 5.05, (0., 2., 0.), no_plot=False)
    # E x B drift: pos_end = pos_ini + time * E/B (if it is a multiple of gyrocycles, which it is.)
    evaluate_test(os.path.join(dir, test_subdir, 'drift_ey.h5'), 6.620796, (0., 1., 0.), no_plot=False)

if __name__ == "__main__":
    run_test()

