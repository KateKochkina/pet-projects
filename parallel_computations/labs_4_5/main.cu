%%writefile main1.cu
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <chrono>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

using namespace std;

int N = 10000;
int BS = 1024;
__device__ const double G = 6.67e-11;
__device__ const double EPS = 1e-15;


template<typename T>
ostream& operator<<(ostream& stream, const vector<T>& vec) {
    for (auto it = vec.begin(); it != vec.end(); ++it) {
        stream << " " << *it;
    }
    return stream;
}

__global__ void plus_product(double* result, double* lhs, double tau, double* rhs, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        result[idx] = lhs[idx] + tau * rhs[idx];
    }
}

__global__ void plus_equal_product(double* result, double tau, double* rhs, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        result[idx] += tau * rhs[idx];
    }
}

//__global__ void calc_k(double* result, double* pos_vel, double* masses, int n) {
//    int idx = blockIdx.x * blockDim.x + threadIdx.x;
//    if (idx < n) {
//        result[3 * idx + 0] = pos_vel[3 * n + 3 * idx + 0];
//        result[3 * idx + 1] = pos_vel[3 * n + 3 * idx + 1];
//        result[3 * idx + 2] = pos_vel[3 * n + 3 * idx + 2];
//        double a_x = 0.0, a_y = 0.0, a_z = 0.0;
//        double dist3, coef;
//        for (int j = 0; j < n; ++j) {
//            dist3 = pow(pow(pos_vel[3 * idx + 0] - pos_vel[3 * j + 0], 2) +
//                                pow(pos_vel[3 * idx + 1] - pos_vel[3 * j + 1], 2) +
//                                pow(pos_vel[3 * idx + 2] - pos_vel[3 * j + 2], 2),
//                        3.0 / 2);
//            coef = masses[j] / max(dist3, EPS);
//            a_x += coef * (pos_vel[3 * idx + 0] - pos_vel[3 * j + 0]);
//            a_y += coef * (pos_vel[3 * idx + 1] - pos_vel[3 * j + 1]);
//            a_z += coef * (pos_vel[3 * idx + 2] - pos_vel[3 * j + 2]);
//        }
//        result[3 * n + 3 * idx + 0] = -G * a_x;
//        result[3 * n + 3 * idx + 1] = -G * a_y;
//        result[3 * n + 3 * idx + 2] = -G * a_z;
//    }
//}

__device__ double sqr(double x) {
    return x * x;
}

__device__ double cube(double x) {
    return x * x * x;
}

__global__ void calc_k(double* result, double* pos_vel, double* masses, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        int n3 = 3 * n;
        int idx3 = 3 * idx;
        result[idx3] = pos_vel[n3 + idx3];
        result[idx3 + 1] = pos_vel[n3 + idx3 + 1];
        result[idx3 + 2] = pos_vel[n3 + idx3 + 2];

        /* copy all positions to shared memory */
        extern __shared__ double pos[];
        pos[idx3] = pos_vel[idx3];
        pos[idx3 + 1] = pos_vel[idx3 + 1];
        pos[idx3 + 2] = pos_vel[idx3 + 2];
        /* copy all masses to shared memory */
        extern __shared__ double mas[];
        mas[idx] = masses[idx];
        __syncthreads();
        /* copy [idx] position to register memory */
        double pos_idx[3] = {pos[idx3], pos[idx3 + 1], pos[idx3 + 2]};
        /* initialize [idx] acceleration on register memory */
        double acc[3] = {0.0, 0.0, 0.0};
        /* initialize coefs on register memory */
        double dist3, coef, diff[3];

        for (int j = 0; j < n; ++j) {
            diff[0] = pos[3 * j] - pos_idx[0];
            diff[1] = pos[3 * j + 1] - pos_idx[1];
            diff[2] = pos[3 * j + 2] - pos_idx[2];
            dist3 = cube(sqrt(sqr(diff[0]) + sqr(diff[1]) + sqr(diff[2])));
            coef = mas[j] / max(dist3, EPS);
            acc[0] += coef * diff[0];
            acc[1] += coef * diff[1];
            acc[2] += coef * diff[2];
        }
        result[n3 + idx3 + 0] = G * acc[0];
        result[n3 + idx3 + 1] = G * acc[1];
        result[n3 + idx3 + 2] = G * acc[2];
    }
}

void runge_kutta_2(const vector<double>& init, const vector<double>& masses, double t_max, double tau,
                   int print_every_n, const string& file_name) {
    ofstream stream(file_name);
    int n = masses.size();
    int size = init.size();  // 6 * n
    int step_num = 0;
    double t_cur = 0.0;
    vector<double> p_cur(init);
    stream << setprecision(10) << t_cur << " " << p_cur << '\n';

    double *d_p_cur, *d_p_mid, *d_masses, *d_k;
    cudaMalloc(&d_p_cur, size * sizeof(double));
    cudaMalloc(&d_p_mid, size * sizeof(double));
    cudaMalloc(&d_masses, n * sizeof(double));
    cudaMalloc(&d_k, size * sizeof(double));

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    float elapsed_time;
    float exec_time = 0.0, copy_time = 0.0;

    cudaEventRecord(start);
    /* copy masses and p(t) to device */
    cudaMemcpy(d_p_cur, p_cur.data(), size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_masses, masses.data(), n * sizeof(double), cudaMemcpyHostToDevice);
    /* end of record */
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsed_time, start, stop);
    copy_time += elapsed_time;

    while (t_cur < t_max) {
        t_cur += tau;

        cudaEventRecord(start);
        /* k(t) = { v(t), a(t) } */
        calc_k<<<(n + BS) / BS, BS>>>(d_k, d_p_cur, d_masses, n);
        /* p(t + tau/2) = p(t) + tau / 2 * k(t) */
        plus_product<<<(size + BS) / BS, BS>>>(d_p_mid, d_p_cur, tau / 2, d_k, size);
        /* k(t + tau/2) = { v(t + tau/2), a(t + tau/2) } */
        calc_k<<<(n + BS) / BS, BS>>>(d_k, d_p_mid, d_masses, n);
        /* p(t + tau) += tau * k(t + tau/2) */
        plus_equal_product<<<(size + BS) / BS, BS>>>(d_p_cur, tau, d_k, size);
        /* end of record */
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&elapsed_time, start, stop);
        exec_time += elapsed_time;

        cudaEventRecord(start);
        /* copy p(t + tau) to host */
        cudaMemcpy(p_cur.data(), d_p_cur, size * sizeof(double), cudaMemcpyDeviceToHost);
        /* end of record */
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&elapsed_time, start, stop);
        copy_time += elapsed_time;

        if (++step_num % print_every_n == 0) {
            stream << setprecision(10) << t_cur << " " << p_cur << '\n';
        }
    }

    cudaDeviceSynchronize();

    cudaFree(d_p_cur);
    cudaFree(d_p_mid);
    cudaFree(d_masses);
    cudaFree(d_k);
//    cout << "Time exec: " << exec_time / 1000.0 << endl;
//    cout << "Time copy: " << copy_time / 1000.0 << endl;
}

int main(int argc, char** argv) {
    // $./run N bs
    if (argc > 1) {
        N = stoi(argv[1]);
    }
    if (argc > 2) {
        BS = stoi(argv[2]);
    }
    cout << "N: " << N << "\nbs: " << BS << "\n";

    /* part 1 data init */
    vector<double> input_data = {1.0, 0.0, 0.0, 0.0, 0.9, 0.0,
                                 0.0, 1.0, 0.0, -0.9, 0.0, 0.0,
                                 -1.0, 0.0, 0.0, 0.0, -0.9, 0.0,
                                 0.0, -1.0, 0.0, 0.9, 0.0, 0.0};
    vector<double> masses = {8810324116.227, 8810324116.227, 8810324116.227, 8810324116.227};
    vector<double> init(input_data.size());
    for (int i = 0; i < init.size() / 6; ++i) {
        init[3 * i + 0] = input_data[6 * i + 0];
        init[3 * i + 1] = input_data[6 * i + 1];
        init[3 * i + 2] = input_data[6 * i + 2];
        init[init.size() / 2 + 3 * i + 0] = input_data[6 * i + 3];
        init[init.size() / 2 + 3 * i + 1] = input_data[6 * i + 4];
        init[init.size() / 2 + 3 * i + 2] = input_data[6 * i + 5];
    }

    /* part 1 calc */
    runge_kutta_2(init, masses, 20, 0.01, 10, "./cuda.txt");

    /* distributions init */
    std::random_device seeder;
    const auto seed = seeder.entropy() ? seeder() : time(nullptr);
    std::mt19937 rnd_gen(static_cast<std::mt19937::result_type>(seed));
    std::uniform_real_distribution<double> distr_pos(-10.0, 10.0);
    std::uniform_real_distribution<double> distr_vel(-1.0, 1.0);
    std::uniform_real_distribution<double> distr_mass(0.1, 10.0);
    auto rnd_pos = std::bind(distr_pos, rnd_gen);
    auto rnd_vel = std::bind(distr_vel, rnd_gen);
    auto rnd_mas = std::bind(distr_mass, rnd_gen);

    /* part 2 data init */
    init.resize(6 * N);
    for (int i = 0; i < 3 * N; ++i) {
        init[i] = rnd_pos();
    }
    for (int i = 3 * N; i < 6 * N; ++i) {
        init[i] = rnd_vel();
    }
    masses.resize(N);
    for (int i = 0; i < N; ++i) {
        masses[i] = rnd_mas();
    }

    /* part 2 calc */
    auto t1 = chrono::high_resolution_clock::now();
    runge_kutta_2(init, masses, 0.1, 0.05, 1, "./cuda_10000.txt");
    auto t2 = chrono::high_resolution_clock::now();
    cout << "Time cuda: " << chrono::duration_cast<chrono::milliseconds>(t2 - t1).count() / 1000.0 << endl;

    return 0;
}