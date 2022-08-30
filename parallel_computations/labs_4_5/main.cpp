#include "mpi.h"
#include "omp.h"
#include <chrono>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

using namespace std;

const int N = 10000;
const int NUM_THREADS = omp_get_num_threads();
const double G = 6.67e-11;

MPI_Datatype MPI_COORDINATES;
MPI_Datatype MPI_POINT;
MPI_Datatype MPI_POINT_POS, _MPI_POINT_POS_;


struct Coordinates {
    double data[3];

    Coordinates()
        : data{.0, .0, .0} {}

    Coordinates(double x1, double x2, double x3)
        : data{x1, x2, x3} {}

    Coordinates(const Coordinates& coord)
        : data{coord.data[0], coord.data[1], coord.data[2]} {}

    friend ostream& operator<<(ostream& stream, const Coordinates& coord);
};

Coordinates operator+(const Coordinates& lhs, const Coordinates& rhs) {
    return Coordinates(lhs.data[0] + rhs.data[0],
                       lhs.data[1] + rhs.data[1],
                       lhs.data[2] + rhs.data[2]);
}

Coordinates operator-(const Coordinates& lhs, const Coordinates& rhs) {
    return Coordinates(lhs.data[0] - rhs.data[0],
                       lhs.data[1] - rhs.data[1],
                       lhs.data[2] - rhs.data[2]);
}

Coordinates operator*(double lhs, const Coordinates& rhs) {
    return Coordinates(lhs * rhs.data[0],
                       lhs * rhs.data[1],
                       lhs * rhs.data[2]);
}

Coordinates& operator+=(Coordinates& lhs, const Coordinates& rhs) {
    for (int i = 0; i < 3; ++i) {
        lhs.data[i] += rhs.data[i];
    }
    return lhs;
}

ostream& operator<<(ostream& stream, const Coordinates& coord) {
    for (int i = 0; i < 3; ++i) {
        stream << coord.data[i] << " ";
    }
    return stream;
}

double get_dist3(const Coordinates& c1, const Coordinates& c2) {
    return pow(pow(c1.data[0] - c2.data[0], 2) +
                       pow(c1.data[1] - c2.data[1], 2) +
                       pow(c1.data[2] - c2.data[2], 2),
               3.0 / 2);
}


struct Point {
    Coordinates x;
    Coordinates v;
    double mass;

    Point()
        : mass(.0) {}

    Point(double x1, double x2, double x3, double v1, double v2, double v3, double mass)
        : x(x1, x2, x3), v(v1, v2, v3), mass(mass) {}

    Point(const Point& point)
        : x(point.x), v(point.v), mass(point.mass) {}
};

ostream& operator<<(ostream& stream, const vector<Point>& vec) {
    for (auto it = vec.begin(); it != vec.end(); ++it) {
        stream << it->x;
    }
    for (auto it = vec.begin(); it != vec.end(); ++it) {
        stream << it->v;
    }
    return stream;
}


vector<Coordinates> calc_accelerations(const vector<Point>& pts, const vector<int>& counts, const vector<int>& displs,
                                       int rank, int np, int nth) {
    double eps = 1e-15;
    size_t n = pts.size();
    vector<Coordinates> accels(counts[rank]);
    vector<Coordinates> accels_global(n);
#pragma omp parallel for num_threads(nth)
    for (int i = displs[rank]; i < displs[rank + 1]; ++i) {
        Coordinates accel(.0, .0, .0);
        for (int j = 0; j < n; ++j) {
            double dist3 = max(get_dist3(pts[i].x, pts[j].x), eps);
            accel += (pts[j].mass / dist3) * (pts[j].x - pts[i].x);
        }
        accels[i - displs[rank]] = G * accel;
    }
    MPI_Gatherv(accels.data(), counts[rank], MPI_COORDINATES,
                accels_global.data(), counts.data(), displs.data(), MPI_COORDINATES, 0, MPI_COMM_WORLD);
    return accels_global;
}

void runge_kutta_2(const vector<Point>& init, double t_max, double tau, int print_every_n,
                   const string& file_name, int rank, int np, int nth) {
    ofstream stream(file_name);
    int n = init.size();
    int step_num = 0;
    double t_cur = 0;
    vector<Point> p_cur(init), p_mid(init);
    if (rank == 0) {
        stream << setprecision(10) << t_cur << " " << p_cur << '\n';
    }
    vector<int> counts(np, n / np);
    for (int i = 0; i < n % np; ++i) {
        ++counts[i];
    }
    vector<int> displs(np + 1, 0);
    for (int i = 0; i < np; ++i) {
        displs[i + 1] = displs[i] + counts[i];
    }
    while (t_cur < t_max) {
        t_cur += tau;
        const vector<Coordinates>& accels1 = calc_accelerations(p_cur, counts, displs, rank, np, nth);
        if (rank == 0) {
            for (int i = 0; i < n; ++i) {
                p_mid[i].x = p_cur[i].x + tau / 2 * p_cur[i].v;
                p_mid[i].v = p_cur[i].v + tau / 2 * accels1[i];
            }
        }
        MPI_Bcast(p_mid.data(), n, MPI_POINT_POS, 0, MPI_COMM_WORLD);

        const vector<Coordinates>& accels2 = calc_accelerations(p_mid, counts, displs, rank, np, nth);
        if (rank == 0) {
            for (int i = 0; i < n; ++i) {
                p_cur[i].x += tau * p_mid[i].v;
                p_cur[i].v += tau * accels2[i];
            }
            if (++step_num % print_every_n == 0) {
                stream << setprecision(10) << t_cur << " " << p_cur << '\n';
            }
        }
        MPI_Bcast(p_cur.data(), n, MPI_POINT_POS, 0, MPI_COMM_WORLD);
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int np, rank, nth;
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    nth = NUM_THREADS;
    if (rank == 0) {
        cout << "np = " << np << "\nnth = " << nth << "\nN = " << N << "\n\n";
    }

    /* create coordinates type */
    const int c_num_fields = 1;
    int c_lengths[c_num_fields] = {3};
    MPI_Aint c_displacements[c_num_fields] = {offsetof(Coordinates, data)};
    MPI_Datatype c_types[c_num_fields] = {MPI_DOUBLE};
    MPI_Type_create_struct(c_num_fields, c_lengths, c_displacements, c_types, &MPI_COORDINATES);
    MPI_Type_commit(&MPI_COORDINATES);

    /* create point type */
    const int p_num_fields = 3;
    int p_lengths[p_num_fields] = {1, 1, 1};
    MPI_Aint p_displacements[p_num_fields] = {offsetof(Point, x), offsetof(Point, v), offsetof(Point, mass)};
    MPI_Datatype p_types[p_num_fields] = {MPI_COORDINATES, MPI_COORDINATES, MPI_DOUBLE};
    MPI_Type_create_struct(p_num_fields, p_lengths, p_displacements, p_types, &MPI_POINT);
    MPI_Type_commit(&MPI_POINT);

    /* create point position type */
    const int pp_num_fields = 1;
    int pp_lengths[pp_num_fields] = {1};
    MPI_Aint pp_displacements[pp_num_fields] = {offsetof(Point, x)};
    MPI_Datatype pp_types[pp_num_fields] = {MPI_COORDINATES};
    MPI_Type_create_struct(pp_num_fields, pp_lengths, pp_displacements, pp_types, &_MPI_POINT_POS_);
    MPI_Type_commit(&_MPI_POINT_POS_);
    MPI_Type_create_resized(_MPI_POINT_POS_, offsetof(Point, x), sizeof(Point), &MPI_POINT_POS);
    MPI_Type_commit(&MPI_POINT_POS);

    /* part 1 data init */
    vector<double> input_data = {1.0, 0.0, 0.0, 0.0, 0.9, 0.0,
                                 0.0, 1.0, 0.0, -0.9, 0.0, 0.0,
                                 -1.0, 0.0, 0.0, 0.0, -0.9, 0.0,
                                 0.0, -1.0, 0.0, 0.9, 0.0, 0.0};
    vector<double> masses = {8810324116.227, 8810324116.227, 8810324116.227, 8810324116.227};
    vector<Point> points_init;
    for (int i = 0; i < input_data.size() / 6; ++i) {
        points_init.emplace_back(input_data[6 * i + 0], input_data[6 * i + 1], input_data[6 * i + 2],
                                 input_data[6 * i + 3], input_data[6 * i + 4], input_data[6 * i + 5],
                                 masses[i]);
    }

    /* part 1 calc */
    runge_kutta_2(points_init, 20, 0.01, 10,
                  np == 1 ? "./output/seq.txt" : "./output/mpi.txt", rank, np, nth);
    MPI_Barrier(MPI_COMM_WORLD);

    /* distributions init */
    std::random_device seeder;
    const auto seed = seeder.entropy() ? seeder() : time(nullptr);
    std::mt19937 rnd_gen(static_cast<std::mt19937::result_type>(seed));
    std::uniform_real_distribution<double> distr_pos(-10.0, 10.0);
    std::uniform_real_distribution<double> distr_vel(-1.0, 1.0);
    std::uniform_real_distribution<double> distr_mass(0.1, 10.0);
    auto rnd_pos = std::bind(distr_pos, rnd_gen);
    auto rnd_vel = std::bind(distr_vel, rnd_gen);
    auto rnd_mass = std::bind(distr_mass, rnd_gen);

    /* part 2 data init */
    points_init.resize(N);
    if (rank == 0) {
        /*
        for (int i = 0; i < N; ++i) {
            points_init[i] = Point(rnd_pos(), rnd_pos(), rnd_pos(), rnd_vel(), rnd_vel(), rnd_vel(), rnd_mass());
        }
        ofstream file("./input/data_10000.txt");
        if (file.is_open()) {
            for (int i = 0; i < N; ++i) {
                file << points_init[i].x << points_init[i].v << points_init[i].mass << '\n';
            }
            file.close();
        }
        */
        ifstream file("./input/data_10000.txt");
        if (file.is_open()) {
            double x1, x2, x3, v1, v2, v3, mass;
            for (int i = 0; i < N; ++i) {
                file >> x1 >> x2 >> x3 >> v1 >> v2 >> v3 >> mass;
                points_init[i] = Point(x1, x2, x3, v1, v2, v3, mass);
            }
        }
    }
    MPI_Bcast(points_init.data(), points_init.size(), MPI_POINT, 0, MPI_COMM_WORLD);

    /* part 2 calc */
    auto t1 = chrono::high_resolution_clock::now();
    runge_kutta_2(points_init, 0.1, 0.05, 1,
                  np == 1 ? "./output/seq_10000.txt" : "./output/mpi_10000.txt", rank, np, nth);
    auto t2 = chrono::high_resolution_clock::now();
    if (rank == 0) {
        cout << (np == 1 ? "Time seq: " : "Time mpi + omp: ") << chrono::duration_cast<chrono::milliseconds>(t2 - t1).count() / 1000.0 << endl;
    }

    MPI_Type_free(&MPI_POINT_POS);
    MPI_Type_free(&MPI_POINT);
    MPI_Type_free(&MPI_COORDINATES);
    MPI_Finalize();

    return 0;
}