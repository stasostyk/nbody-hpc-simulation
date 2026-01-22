import matplotlib.pyplot as plt
import numpy as np
from get_results_data import *


def extract_step10_data(results_dict, steps):
    n_values = []
    times = []
    for key, time in results_dict.items():
        steps_str = key.split()[1].split('=')[1]
        if int(steps_str) == steps:
            # print(f'key: {key}')
            # Extract n value
            n_str = key.split()[0].split('=')[1]
            n = int(n_str)
            if n <= 100 and steps <= 10:
                continue
            # if n >= int(1e6):
            #     continue
            # Only include non-zero times for better visualization
            # if time > 0:
            n_values.append(n)
            times.append(max(time, 1)) # clamp for better visualization
    # Sort by n
    sorted_pairs = sorted(zip(n_values, times))
    return [p[0] for p in sorted_pairs], [p[1] for p in sorted_pairs]



def plot_results():
    steps = 10

    # Extract data for all algorithms first
    serial_n, serial_time = extract_step10_data(serialResults2D, steps)
    serial_reduced_n, serial_reduced_time = extract_step10_data(serialReducedResults2D, steps)
    serial_omp_n, serial_omp_time = extract_step10_data(serialResults2DOMP4, steps)
    serial_reduced_omp_n, serial_reduced_omp_time = extract_step10_data(serialReducedResults2DOMP4, steps)
    bh_n, bh_time = extract_step10_data(bh2Dresults052D, steps)
    mpi_n, mpi_time = extract_step10_data(mpi2Dresults_n4, steps)
    mpi_reduced_n, mpi_reduced_time = extract_step10_data(mpiReduced2Dresults_n4, steps)
    bh_3D_n, bh_3D_time = extract_step10_data(bh3Dresults05, steps)

    # Create the plot
    plt.figure(figsize=(10, 6))

    # print(f'bh, n: {bh_n}, times: {bh_time}')

    # Plot all algorithms
    # plt.loglog(serial_n, serial_time, 'o-', label='Serial', linewidth=2, markersize=8)
    # plt.loglog(serial_reduced_n, serial_reduced_time, 's-', label='Serial Reduced', linewidth=2, markersize=8)
    # plt.loglog(serial_omp_n, serial_omp_time, 'v-', label='Serial OMP (t=4)', linewidth=2, markersize=8)
    # plt.loglog(serial_reduced_omp_n, serial_reduced_omp_time, '^-', label='Serial Reduced OMP (t=4)', linewidth=2, markersize=8)
    # plt.loglog(mpi_n, mpi_time, 'p-', label='MPI (n=4)', linewidth=2, markersize=8)
    # plt.loglog(mpi_reduced_n, mpi_reduced_time, '*-', label='MPI Reduced (n=4)', linewidth=2, markersize=8)
    plt.loglog(bh_n, bh_time, 'd-', label='2D Barnes-Hut (θ=0.5)', linewidth=2, markersize=8)
    # plt.semilogx(bh_n, bh_time, 'd-', label='2D Barnes-Hut (θ=0.5)', linewidth=2, markersize=8)
    plt.loglog(bh_3D_n, bh_3D_time, 'p-', label='3D Barnes-Hut (θ=0.5)', linewidth=2, markersize=8)
    # plt.semilogx(bh_3D_n, bh_3D_time, 'p-', label='3D Barnes-Hut (θ=0.5)', linewidth=2, markersize=8)

    plt.xlabel('Number of particles (n)', fontsize=12)
    plt.ylabel('Execution time (ms)', fontsize=12)
    plt.title(f'N-body Simulation Performance Comparison (steps={steps})', fontsize=14, fontweight='bold')
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.3, which='both')

    # Make the plot look nicer
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    plot_results()