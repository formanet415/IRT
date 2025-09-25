import numpy as np
import h5py
import matplotlib.pyplot as plt


# TODO: implement a pass/fail test condition based on an analytical formula

def plot_test_particle(filename):
    uniform_bz = h5py.File(filename)  # Modify this for your case

    x  = uniform_bz["x"][:]
    vx = uniform_bz["vx"][:]
    vy = uniform_bz["vy"][:]
    vz = uniform_bz["vz"][:]
    nstep = np.arange(len(x))

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
    plt.show()

# Call the function for each file
plot_test_particle("/Users/formanet/github/MCHPC/build/tests/boris/uniform_bz.h5")
plot_test_particle("/Users/formanet/github/MCHPC/build/tests/boris/drift_ey.h5")

