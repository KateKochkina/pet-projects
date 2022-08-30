#include "mpi.h"
#include <algorithm>
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <set>
#include <vector>

using std::cout;
using std::endl;
using std::set;
using std::vector;

const size_t n = 100;
const double h = 1.0 / n;
const double k = 10.0;
const double EPS = 1e-4;

/* methods
 */
const char JACOBI[] = "Jacobi";
const char ZEIDEL[] = "Zeidel";

/* 1 - MPI_Send, MPI_Recv
 * 2 - MPI_Sendrecv
 * 3 - MPI_Send_init, MPI_Recv_init
 */
int SEND_TYPE_DEFAULT = 1;

/* tags
 */
int SEND_LEFT = 1;
int SEND_RIGHT = 2;


double f(double x, double y) {
    return (2 + (k * k + M_PI * M_PI) * (1 - x) * x) * sin(M_PI * y);
}

double u(double x, double y) {
    return (1 - x) * x * sin(M_PI * y);
}

double norm(double* vec, size_t n) {
    double sum = 0.0;
    for (size_t i = 0; i < n; ++i) {
        sum += fabs(vec[i]);
    }
    return sum;
}

void method(const char name[], int np, int rank, size_t send_type) {
    int y_size = n + 1;
    int x_size;
    vector<double> sol, sol_next, delta, rhs;
    vector<double> sol_global, delta_global, rhs_global, exact_solution;
    vector<int> x_size_vec, sendcounts, displs;

    if (rank == 0) {
        sol_global.resize(y_size * y_size);
        delta_global.resize(y_size * y_size);
        vector<double> points(y_size);
        for (size_t i = 0; i < y_size; ++i) {
            points[i] = i * h;
        }
        rhs_global.resize(y_size * y_size);
        exact_solution.resize(y_size * y_size);
        for (size_t i = 0; i < y_size; ++i) {
            for (size_t j = 0; j < y_size; ++j) {
                rhs_global[i * y_size + j] = h * h * f(points[i], points[j]);
                exact_solution[i * y_size + j] = u(points[i], points[j]);
            }
        }
        x_size_vec.resize(np, y_size / np);
        for (int i = 0; i < y_size % np; ++i) {
            ++x_size_vec[i];
        }
        sendcounts.resize(np);
        for (int i = 0; i < np; ++i) {
            sendcounts[i] = x_size_vec[i] * y_size;
        }
        displs.resize(np);
        displs[0] = 0;
        for (int i = 1; i < np; ++i) {
            displs[i] = (displs[i - 1] / y_size + x_size_vec[i - 1]) * y_size;
        }
    }
    MPI_Scatter(x_size_vec.data(), 1, MPI_INT, &x_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    x_size += 2;  /* add ghost cells */
    sol.resize(x_size * y_size, 0.0);
    sol_next.resize(x_size * y_size, 0.0);
    if (np > 1) {
        delta.resize(x_size * y_size);
    }
    rhs.resize(x_size * y_size);
    MPI_Scatterv(rhs_global.data(), sendcounts.data(), displs.data(), MPI_DOUBLE,
                 rhs.data() + y_size, (x_size - 2) * y_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    set<double> boundary_index;  /* all boundary points indexes */
    if (rank == 0) {
        for (int i = 0; i < y_size; ++i) {
            boundary_index.insert(y_size + i);
        }
        for (int i = 1; i < x_size; ++i) {
            boundary_index.insert(i * y_size);
            boundary_index.insert(i * y_size + y_size - 1);
        }
        if (np == 1) {
            for (int i = 0; i < y_size; ++i) {
                boundary_index.insert((x_size - 2) * y_size + i);
            }
        }
    } else if (rank == np - 1) {
        for (int i = 0; i < y_size; ++i) {
            boundary_index.insert((x_size - 2) * y_size + i);
        }
        for (int i = 0; i < x_size - 1; ++i) {
            boundary_index.insert(i * y_size);
            boundary_index.insert(i * y_size + y_size - 1);
        }
    } else {
        for (int i = 0; i < x_size; ++i) {
            boundary_index.insert(i * y_size);
            boundary_index.insert(i * y_size + y_size - 1);
        }
    }

    double t1 = MPI_Wtime();

    int send_left_cnt = (rank == 0) ? 0 : y_size,
        send_right_cnt = (rank == np - 1) ? 0 : y_size,
        neib_left = (rank == 0) ? np - 1 : rank - 1,
        neib_right = (rank == np - 1) ? 0 : rank + 1;
    MPI_Status st;
    MPI_Status *stSend = new MPI_Status[2],
               *stRecv = new MPI_Status[2];
    MPI_Request *reqSend1 = new MPI_Request[2],
                *reqRecv1 = new MPI_Request[2],
                *reqSend2 = new MPI_Request[2],
                *reqRecv2 = new MPI_Request[2];
    if (np > 1 && send_type == 3) {
        MPI_Send_init(sol.data() + (x_size - 2) * y_size, send_right_cnt, MPI_DOUBLE, neib_right, SEND_RIGHT, MPI_COMM_WORLD, reqSend1);
        MPI_Recv_init(sol.data(), send_left_cnt, MPI_DOUBLE, neib_left, SEND_RIGHT, MPI_COMM_WORLD, reqRecv1);
        MPI_Send_init(sol.data() + y_size, send_left_cnt, MPI_DOUBLE, neib_left, SEND_LEFT, MPI_COMM_WORLD, reqSend1 + 1);
        MPI_Recv_init(sol.data() + (x_size - 1) * y_size, send_right_cnt, MPI_DOUBLE, neib_right, SEND_LEFT, MPI_COMM_WORLD, reqRecv1 + 1);

        MPI_Send_init(sol_next.data() + (x_size - 2) * y_size, send_right_cnt, MPI_DOUBLE, neib_right, SEND_RIGHT, MPI_COMM_WORLD, reqSend2);
        MPI_Recv_init(sol_next.data(), send_left_cnt, MPI_DOUBLE, neib_left, SEND_RIGHT, MPI_COMM_WORLD, reqRecv2);
        MPI_Send_init(sol_next.data() + y_size, send_left_cnt, MPI_DOUBLE, neib_left, SEND_LEFT, MPI_COMM_WORLD, reqSend2 + 1);
        MPI_Recv_init(sol_next.data() + (x_size - 1) * y_size, send_right_cnt, MPI_DOUBLE, neib_right, SEND_LEFT, MPI_COMM_WORLD, reqRecv2 + 1);
    }

    double denominator = 4 + h * h * k * k;
    double sol_norm, delta_norm;
    double sol_norm_global, delta_norm_global;
    double error = 1e10;
    int iter = 0;

    while (error > EPS) {
        ++iter;
        if (np > 1) {
            if (send_type == 1) {
                MPI_Send(sol.data() + (x_size - 2) * y_size, send_right_cnt, MPI_DOUBLE, neib_right, SEND_RIGHT, MPI_COMM_WORLD);
                MPI_Recv(sol.data(), send_left_cnt, MPI_DOUBLE, neib_left, SEND_RIGHT, MPI_COMM_WORLD, &st);
                MPI_Send(sol.data() + y_size, send_left_cnt, MPI_DOUBLE, neib_left, SEND_LEFT, MPI_COMM_WORLD);
                MPI_Recv(sol.data() + (x_size - 1) * y_size, send_right_cnt, MPI_DOUBLE, neib_right, SEND_LEFT, MPI_COMM_WORLD, &st);
            } else if (send_type == 2) {
                MPI_Sendrecv(sol.data() + (x_size - 2) * y_size, send_right_cnt, MPI_DOUBLE, neib_right, SEND_RIGHT,
                             sol.data(), send_left_cnt, MPI_DOUBLE, neib_left, SEND_RIGHT, MPI_COMM_WORLD, &st);
                MPI_Sendrecv(sol.data() + y_size, send_left_cnt, MPI_DOUBLE, neib_left, SEND_LEFT,
                             sol.data() + (x_size - 1) * y_size, send_right_cnt, MPI_DOUBLE, neib_right, SEND_LEFT, MPI_COMM_WORLD, &st);
            } else if (send_type == 3) {
                if (iter % 2 == 1) {
                    MPI_Startall(2, reqSend1);
                    MPI_Startall(2, reqRecv1);
                } else {
                    MPI_Startall(2, reqSend2);
                    MPI_Startall(2, reqRecv2);
                }
            }
        }

        if (name == JACOBI) {
            /* all except left and right bounds */
            for (int i = 2 * y_size; i < (x_size - 2) * y_size; ++i) {
                if (boundary_index.count(i) == 0) {
                    sol_next[i] = (rhs[i] + sol[i - 1] + sol[i + 1] + sol[i + y_size] + sol[i - y_size]) / denominator;
                }
            }

            if (np > 1 && send_type == 3) {
                if (iter % 2 == 1) {
                    MPI_Waitall(2, reqSend1, stSend);
                    MPI_Waitall(2, reqRecv1, stSend);
                } else {
                    MPI_Waitall(2, reqSend2, stSend);
                    MPI_Waitall(2, reqRecv2, stSend);
                }
            }

            /* left bound */
            for (int i = y_size; i < 2 * y_size; ++i) {
                if (boundary_index.count(i) == 0) {
                    sol_next[i] = (rhs[i] + sol[i - 1] + sol[i + 1] + sol[i + y_size] + sol[i - y_size]) / denominator;
                }
            }

            /* right bound */
            for (int i = (x_size - 2) * y_size; i < (x_size - 1) * y_size; ++i) {
                if (boundary_index.count(i) == 0) {
                    sol_next[i] = (rhs[i] + sol[i - 1] + sol[i + 1] + sol[i + y_size] + sol[i - y_size]) / denominator;
                }
            }

        } else if (name == ZEIDEL) {
            if (np > 1 && send_type == 3) {
                if (iter % 2 == 1) {
                    MPI_Waitall(2, reqSend1, stSend);
                    MPI_Waitall(2, reqRecv1, stSend);
                } else {
                    MPI_Waitall(2, reqSend2, stSend);
                    MPI_Waitall(2, reqRecv2, stSend);
                }
            }
            
            /* even layers with old solution */
            for (int i = y_size; i < (x_size - 1) * y_size; i += 2) {
                if (boundary_index.count(i) == 0) {
                    sol_next[i] = (rhs[i] + sol[i - 1] + sol[i + 1] + sol[i + y_size] + sol[i - y_size]) / denominator;
                }
            }

            if (np > 1) {
                if (send_type == 1) {
                    MPI_Send(sol_next.data() + (x_size - 2) * y_size, send_right_cnt, MPI_DOUBLE, neib_right, SEND_RIGHT, MPI_COMM_WORLD);
                    MPI_Recv(sol_next.data(), send_left_cnt, MPI_DOUBLE, neib_left, SEND_RIGHT, MPI_COMM_WORLD, &st);
                    MPI_Send(sol_next.data() + y_size, send_left_cnt, MPI_DOUBLE, neib_left, SEND_LEFT, MPI_COMM_WORLD);
                    MPI_Recv(sol_next.data() + (x_size - 1) * y_size, send_right_cnt, MPI_DOUBLE, neib_right, SEND_LEFT, MPI_COMM_WORLD, &st);
                } else if (send_type == 2) {
                    MPI_Sendrecv(sol_next.data() + (x_size - 2) * y_size, send_right_cnt, MPI_DOUBLE, neib_right, SEND_RIGHT,
                                 sol_next.data(), send_left_cnt, MPI_DOUBLE, neib_left, SEND_RIGHT, MPI_COMM_WORLD, &st);
                    MPI_Sendrecv(sol_next.data() + y_size, send_left_cnt, MPI_DOUBLE, neib_left, SEND_LEFT,
                                 sol_next.data() + (x_size - 1) * y_size, send_right_cnt, MPI_DOUBLE, neib_right, SEND_LEFT, MPI_COMM_WORLD, &st);
                } else if (send_type == 3) {
                    if (iter % 2 == 1) {
                        MPI_Startall(2, reqSend2);
                        MPI_Startall(2, reqRecv2);
                        MPI_Waitall(2, reqSend2, stSend);
                        MPI_Waitall(2, reqRecv2, stSend);
                    } else {
                        MPI_Startall(2, reqSend1);
                        MPI_Startall(2, reqRecv1);
                        MPI_Waitall(2, reqSend1, stSend);
                        MPI_Waitall(2, reqRecv1, stSend);
                    }
                }
            }

            /* odd layers with new solution */
            for (int i = y_size + 1; i < (x_size - 1) * y_size; i += 2) {
                if (boundary_index.count(i) == 0) {
                    sol_next[i] = (rhs[i] + sol_next[i - 1] + sol_next[i + 1] + sol_next[i + y_size] + sol_next[i - y_size]) / denominator;
                }
            }
        }

        sol.swap(sol_next);

        if (np > 1) {
            for (int i = 0; i < x_size * y_size; ++i) {
                delta[i] = fabs(sol_next[i] - sol[i]);
            }
            sol_norm = norm(sol.data(), sol.size());
            MPI_Reduce(&sol_norm, &sol_norm_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            delta_norm = norm(delta.data(), delta.size());
            MPI_Reduce(&delta_norm, &delta_norm_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        } else {
            for (int i = 0; i < y_size * y_size; ++i) {
                sol_global[i] = sol[y_size + i];
                delta_global[i] = fabs(sol_next[y_size + i] - sol[y_size + i]);
            }
            sol_norm_global = norm(sol_global.data(), sol_global.size());
            delta_norm_global = norm(delta_global.data(), delta_global.size());
        }

        if (rank == 0) {
            error = delta_norm_global / sol_norm_global;
        }
        if (np > 1) {
            MPI_Bcast(&error, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
    }

    if (np > 1) {
        MPI_Gatherv(sol.data() + y_size, (x_size - 2) * y_size, MPI_DOUBLE,
                    sol_global.data(), sendcounts.data(), displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    double t2 = MPI_Wtime();

    if (rank == 0) {
        if (send_type == 1) {
            cout << name << " (MPI_Send + MPI_Recv)\n";
        }
        else if (send_type == 2) {
            cout << name << " (MPI_Sendrecv)\n";
        }
        else if (send_type == 3) {
            cout << name << " (MPI_Send_init + MPI_Recv_init)\n";
        }
        cout << "Number of iterations: " << iter << endl;
        cout << "Time: " << t2 - t1 << endl;
        for (int i = 0; i < y_size * y_size; ++i) {
            delta_global[i] = fabs(sol_global[i] - exact_solution[i]);
        }
        sol_norm_global = norm(exact_solution.data(), exact_solution.size());
        delta_norm_global = norm(delta_global.data(), delta_global.size());
        error = delta_norm_global / sol_norm_global;
        cout << "Error: " << error << endl << endl;
    }
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int np, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        cout << "np: " << np << endl << endl;
    }

    int send_type = SEND_TYPE_DEFAULT;
    /* if (send_type == 1 || send_type == 2 || send_type == 3) */
    for (send_type = 1; send_type <= 3; ++send_type) {
        method(JACOBI, np, rank, send_type);
        method(ZEIDEL, np, rank, send_type);
    }

    MPI_Finalize();
    return 0;
}
