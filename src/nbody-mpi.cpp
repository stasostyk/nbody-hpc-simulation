#include <mpi.h>
#include <assert.h>
#include <vector>
#include <cmath>
#include "adaptive_timestep.hpp"

#include "acceleration-accumulator.hpp"
#include "body.hpp"
#include "forces/gravity.hpp"
#include "integrators/euler.hpp"
#include "integrators/integrator.hpp"
#include "integrators/sympletic.hpp"
#include "integrators/verlet.hpp"
#include "integrators/rk4.hpp"
#include "mpi-accumulator.hpp"
#include "utils.hpp"
<<<<<<< Updated upstream
#include <algorithm>
=======
#include "timer.hpp"
#include "adaptive_dt.hpp"

#if USE_OPENMP
    #include <omp.h>
#endif
>>>>>>> Stashed changes

constexpr int DIM = 3;

static MPI_Datatype MPI_VEC;

inline static void initMPIType() {
    MPI_Type_contiguous(DIM, MPI_DOUBLE, &MPI_VEC);
    MPI_Type_commit(&MPI_VEC);
}

inline static void freeMPIType() {
    MPI_Type_free(&MPI_VEC);
}

// utility: block counts/displacements for an array of N elements across comm_sz
void make_counts_displs(int N, int comm_sz, std::vector<int> &counts, std::vector<int> &displs) {
    counts.resize(comm_sz);
    displs.resize(comm_sz);
    int base = N / comm_sz;
    int rem = N % comm_sz;
    int off = 0;
    for (int i = 0; i < comm_sz; ++i) {
        counts[i] = base + (i < rem ? 1 : 0);
        displs[i] = off;
        off += counts[i];
    }
}

inline void allMPIInit(int *argc, char ***argv, int &mpiSize, int &mpiRank) {
    MPI_Init(argc, argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    initMPIType(); // to be able to use MPI_VEC for Vecs
}

inline void allMPIFinalize() {
    freeMPIType();
    MPI_Finalize();
}

int main(int argc, char** argv) {
    constexpr int outputStride = 1;

    const forces::force<DIM, EmptyAttributes> &force = forces::gravity<DIM>();
    bodies<DIM, EmptyAttributes> bodies;

    int mpiSize, mpiRank;
    allMPIInit(&argc, &argv, mpiSize, mpiRank);
    
    int n;
    int steps;
    double dt;

    if (mpiRank == 0) {
        // Root process reads input, and will broadcast the data.
        utils::readFromFile("test1.in.out", steps, dt, bodies);
        n = bodies.globalSize();
    } 

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&steps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // TODO make this computation only in root process (mpiRank 0)
    std::vector<int> counts, displs;
    make_counts_displs(n, mpiSize, counts, displs);
    int locN = counts[mpiRank];
    int locOffset = displs[mpiRank];

    bodies.resize(n, locN, locOffset);

    MPI_Bcast(bodies.position.data(), n, MPI_VEC, 0, MPI_COMM_WORLD);
    MPI_Bcast(bodies.velocity.data(), n, MPI_VEC, 0, MPI_COMM_WORLD);
    MPI_Bcast(bodies.mass.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    
    MPIAccumulator<DIM, EmptyAttributes> accumulator(MPI_VEC, locN, counts, displs, force);
    // integrators::Euler<DIM, EmptyAttributes> integrator(accumulator);
    // integrators::Sympletic<DIM, EmptyAttributes> integrator(accumulator);
    // integrators::Verlet<DIM, EmptyAttributes> integrator(accumulator);
    integrators::RK4<DIM, EmptyAttributes> integrator(accumulator);
    
    const double dt_max = dt;
    const double T_end  = steps * dt_max;

<<<<<<< Updated upstream
    double time = 0.0;
    double next_output_time = dt_max;
    int frame = 0;

    std::vector<double> particle_dt(bodies.localSize(), dt_max);

    while (time + 1e-15 < T_end) {
        accumulator.compute(bodies);

        double local_dt =
            timestep::update_timesteps<DIM, EmptyAttributes>(
                bodies, accumulator, dt_max, particle_dt);

        double dt_global = local_dt;
        MPI_Allreduce(&local_dt, &dt_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

        // propose adaptive step
        double dt_step = std::min(dt_global, T_end - time);

        // if we'd cross the next output time, land exactly on it
        dt_step = std::min(dt_step, next_output_time - time);

        // safety (floating point edge)
        if (dt_step < 0.0) dt_step = 0.0;

        // advance (only if positive)
        if (dt_step > 0.0) {
            integrator.step(bodies, dt_step);
            time += dt_step;
        } else {
            // already at an output time within tolerance
            time = next_output_time;
        }

        // catch-up outputs
        while (frame < steps && time >= next_output_time - 1e-12 * dt_max) {

            MPI_Allgatherv(
                bodies.position.data() + bodies.localOffset(),
                bodies.localSize(), MPI_VEC,
                bodies.position.data(),
                counts.data(), displs.data(),
                MPI_VEC, MPI_COMM_WORLD
            );

            if (mpiRank == 0) {
                utils::saveToFile("test1-MPI." + std::to_string(frame) + ".out",
                                steps, dt_max, bodies, false);
            }

            next_output_time = std::min(next_output_time + dt_max, T_end);
            frame++;
=======
    const double dt0 = dt;
    const double t_end = dt0 * static_cast<double>(steps);

    const double dt_max = dt0 * 4.0;
    const double dt_min = dt0 * 1e-6;
    const double eps_v = 5e-7;
    const double v_floor = 1e-6;
    const double a_floor = 1e-12;
    const double max_growth = 1.5;

    double time = 0.0;
    const int out_every = 10;
    double next_out = dt0 * out_every;
    int frame = 0;

    accumulator.compute(bodies);

    double dt_prev = dt0;

    const double dt_raw0_local = compute_local_dt<DIM>(bodies, accumulator, dt_prev, 1e100, dt_min, eps_v, v_floor, a_floor, 1e100);

    double dt_raw0 = dt_raw0_local;
    MPI_Allreduce(&dt_raw0_local, &dt_raw0, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    double dt_local = compute_local_dt<DIM>(bodies, accumulator, dt_prev, dt_max, dt_min, eps_v, v_floor, a_floor, max_growth);

    double dt_curr = dt_local;
    MPI_Allreduce(&dt_local, &dt_curr, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    if (mpiRank == 0) {
    std::cout << "time=" << time << " dt_raw=" << dt_raw0 << " dt_curr=" << dt_curr << " dt_max=" << dt_max << std::endl;
    }


    while (time < t_end && !time_reached(time, t_end)) {
    double dt_step = dt_curr;
    dt_step = std::min(dt_step, t_end - time);
    dt_step = std::min(dt_step, next_out - time);

    integrator.step(bodies, dt_step);
    time += dt_step;

    accumulator.compute(bodies);

    double dt_local_next =
        compute_local_dt<DIM>(bodies, accumulator, dt_step, dt_max, dt_min, eps_v, v_floor, a_floor, max_growth);

    double dt_curr_next = dt_local_next;
    MPI_Allreduce(&dt_local_next, &dt_curr_next, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    if (mpiRank == 0) {
        std::cout << "time=" << time
                << " dt_step=" << dt_step
                << " dt_curr=" << dt_curr_next
                << " dt_local(min)=" << dt_local_next
                << std::endl;
    }

    dt_local = dt_local_next;
    dt_curr = dt_curr_next;

    if (time_reached(time, next_out)) {
        if (frame % outputStride == 0) {
        MPI_Allgatherv(bodies.position.data() + bodies.localOffset(),
                        bodies.localSize(), MPI_VEC, bodies.position.data(),
                        counts.data(), displs.data(), MPI_VEC, MPI_COMM_WORLD);
        MPI_Allgatherv(bodies.velocity.data() + bodies.localOffset(),
                        bodies.localSize(), MPI_VEC, bodies.velocity.data(),
                        counts.data(), displs.data(), MPI_VEC, MPI_COMM_WORLD);
        if (mpiRank == 0) {
            utils::saveToFile("test1-MPI." + std::to_string(frame) + ".out",
                            steps, dt0, bodies, false);
>>>>>>> Stashed changes
        }
        }

        ++frame;
        next_out = std::min(t_end, dt0 * static_cast<double>((frame + 1) * out_every));
    }
    }




    if (mpiRank == 0 && frame != steps) {
        std::cerr << "Warning: produced " << frame
                << " frames, expected " << steps << "\n";
    }





    allMPIFinalize();
    return 0;
}
