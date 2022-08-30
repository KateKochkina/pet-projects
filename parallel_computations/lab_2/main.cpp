#include "omp.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <random>

size_t NUM_THREADS = omp_get_max_threads();
size_t N = 3e6;  // number of particles
size_t M = 5e5;  // number of mass calc points
double EPS = 1e-6;


struct Body {
    double x;
    double m;
};

struct Node {
    Body body;
    Node* left;
    Node* right;

    Node(Body body)
        : body(body), left(nullptr), right(nullptr) {}
};

class Tree {
public:
    Tree(Body* body, size_t n) : root(new Node(body[0])) {
#pragma omp parallel for num_threads(NUM_THREADS)
        for (size_t i = 1; i < n; ++i) {
            Node* current_node = root;
            bool is_new_node_added = false;
            while (!is_new_node_added) {
                if (body[i].x < current_node->body.x) {
                    if (current_node->left != nullptr) {
                        current_node = current_node->left;
                    }
                    if (current_node->left == nullptr) {
#pragma omp critical
                        if (current_node->left == nullptr) {
                            is_new_node_added = true;
                            current_node->left = new Node(body[i]);
                        }
                    }
                } else if (body[i].x > current_node->body.x) {
                    if (current_node->right != nullptr) {
                        current_node = current_node->right;
                    }
                    if (current_node->right == nullptr) {
#pragma omp critical
                        if (current_node->right == nullptr) {
                            is_new_node_added = true;
                            current_node->right = new Node(body[i]);
                        }
                    }
                } else {
                    is_new_node_added = true;
                }
            }
        }
    }

    double get_mass_sum(double x, double eps) {
        return this->get_mass_sum_in_subtree(this->root, x, eps);
    }

    ~Tree() {
        this->delete_subtree(this->root);
    }

private:
    Node* root;

    double get_mass_sum_in_subtree(Node* node, double x, double eps) {
        if (node == nullptr) {
            return 0.0;
        }
        /* smaller than left bound */
        if (node->body.x < x - eps) {
            return get_mass_sum_in_subtree(node->right, x, eps);
        }
        /* greater than right bound */
        if (node->body.x > x + eps) {
            return get_mass_sum_in_subtree(node->left, x, eps);
        }
        /* inside an interval */
        double self_mass = 0.0;
        if (fabs(node->body.x - x) > 1e-16) {
            self_mass = node->body.m;
        }
        return get_mass_sum_in_subtree(node->left, x, eps) +
               get_mass_sum_in_subtree(node->right, x, eps) + self_mass;
    }

    void delete_subtree(Node* node) {
        if (node != nullptr) {
            delete_subtree(node->left);
            delete_subtree(node->right);
            delete node;
        }
    }
};

bool is_close(double x, double y) {
    return fabs(x - y) > 1e-16 && fabs(x - y) < EPS;
}

int main() {
    for (size_t nth = 1; nth <= 4; ++nth) {
        NUM_THREADS = nth;
        std::cout << "\nnum threads: " << NUM_THREADS << "\n";
        std::srand(42);

        Body* points = new Body[N];
        Body* mass_points = new Body[M];
        Tree* tree;

        auto begin = std::chrono::steady_clock::now();
        auto end = std::chrono::steady_clock::now();

        for (size_t i = 0; i < N; ++i) {
            points[i] = {
                1000.0 * (2.0 * std::rand() / RAND_MAX - 1.0),
                (double) std::rand() / RAND_MAX
            };
        }
        std::shuffle(points, points + N, std::default_random_engine(42));

        for (size_t i = 0; i < M; ++i)
        {
            mass_points[i] = points[i * (N / M)];
        }

        begin = std::chrono::steady_clock::now();

        tree = new Tree(points, N);

        end = std::chrono::steady_clock::now();
        auto build_tree_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
        std::cout << "build tree time: " << build_tree_time.count() << " ms\n";

        begin = std::chrono::steady_clock::now();

        double summ = 0.0;
#pragma omp parallel for num_threads(NUM_THREADS) reduction(+: summ)
        for (size_t i = 0; i < M; ++i) {
            summ += tree->get_mass_sum(mass_points[i].x, EPS);
        }

        end = std::chrono::steady_clock::now();
        auto sum_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
        std::cout << "sum time: " << sum_time.count() << " ms\n";

        std::cout << "overall time: " << build_tree_time.count() + sum_time.count() << " ms\n";
        std::cout << "sum: " << summ << "\n\n";

        /*
        double sum_test = 0.0;
#pragma omp parallel for num_threads(NUM_THREADS) reduction(+: sum_test)
        for (size_t i = 0; i < M; ++i) {
           for (size_t j = 0; j < N; ++j) {
               if (is_close(points[j].x, mass_points[i].x)) {
                   sum_test += points[j].m;
               }
           }
        }
        std::cout << "sum test: " << sum_test << "\n";
        */

        delete[] points;
        delete[] mass_points;
        delete tree;
    }
}
