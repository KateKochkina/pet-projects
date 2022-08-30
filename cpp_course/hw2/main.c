#include <stdio.h>

#include "matrix.h"

int main(void) {
    Matrix *m, *m1, *m2, *m3, *m4, *m5, *m6, *m7;
    m = create_matrix_from_file("file.txt");
    m1 = create_matrix(3, 3);
    double d;
    size_t s;
    get_rows(m, &s);
    get_cols(m, &s);
    get_elem(m, 3, 4, &d);
    set_elem(m, 3, 2, 2.56);
    m = mul_scalar(m, 1.1);
    view_matrix(m);
    m2 = transp(m1);
    m3 = sum(m1, m2);
    m4 = sub(m1, m3);
    m5 = mul(m1, m4);
    det(m1, &d);
    m6 = adj(m5);
    m7 = inv(m6);
    free_matrix(m);
    free_matrix(m2);
    free_matrix(m3);
    free_matrix(m4);
    free_matrix(m5);
    free_matrix(m6);
    free_matrix(m7);
    return 0;
}
