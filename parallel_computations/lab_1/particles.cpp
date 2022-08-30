#include "stdlib.h"
#include "stdio.h"
#include "mpi.h"
#include <vector>

using std::vector;

uint32_t N_PARTICLES = 1e6;
uint32_t N_POINTS = 1e2;
double EPS = 0.01;
uint32_t N_THREADS = 4;


struct Particle
{
    double x; // положение (на прямой), совпадающих нет
    double m; // масса

}

class FenwickTree {
public:
    FenwickTree(vector<int64_t> arr, uint32_t size)
        : arr(std::move(arr)), size(size) {
        this->tree = vector<int64_t>(size);
        for (int i = 0; i < size; ++i) {
            int64_t sum = 0.0;
            for (int j = F(i); j <= i; ++j) {
                sum += this->arr[j].m;
            }
            this->tree[i] = sum;
        }
    }

    int64_t rsq(uint32_t l, uint32_t r) {
        if (l == 0) {
            return this->sum(r);
        }
        return this->sum(r) - this->sum(l - 1);
    }

    void add(uint32_t i, int64_t x) {
        while (i < this->size) {
            this->tree[i] += x;
            i |= i + 1;
        }
    }

    void set(uint32_t i, int64_t x) {
        int64_t diff = x - this->arr[i];
        this->arr[i] = x;
        this->add(i, diff);
    }

private:
    static uint32_t F(uint32_t i) {
        return i & (i + 1);
    }

    int64_t sum(uint32_t i) {
        int64_t res = 0.0;
        int32_t j = i;
        while (j >= 0) {
            res += this->tree[j];
            j = F(j) - 1;
        }
        return res;
    }

    const uint32_t size;
    vector<int64_t> arr;
    vector<int64_t> tree;
};

struct Tree
{
    Tree *left, *right; // указатели на потомков (nullptr, если ячейка - лист)
    Body pnt; // материальная точка, соответствующая узлу дерева
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc,&argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    vector<int> v;

    printf("size = %d, rank = %d\n", size, rank);

    MPI_Finalize();
    return 0;
}
