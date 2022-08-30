#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "matrix.h"

void view_matrix(const Matrix* matrix) {
    for (int i = 0; i < matrix->rows; ++i) {
        for (int j = 0; j < matrix->cols; ++j) {
            printf("%f ", matrix->matr[i][j]);
        }
        printf("\n");
    }
}

Matrix* create_matrix_from_file(const char* path_file) {
    FILE* f = fopen(path_file, "r");
    if (!f) {
        return NULL;
    }

    int rows, cols;
    int fields_read = fscanf(f, "%d %d", &rows, &cols);
    if (fields_read != 2) {
        fclose(f);
        fprintf(stderr, "Can't read rows/cols\n");
        return NULL;
    }

    Matrix* matrix = create_matrix(rows, cols);
    if (!matrix) {
        fclose(f);
        return NULL;
    }

    for (int i = 0; i < matrix->rows; ++i) {
        for (int j = 0; j < matrix->cols; ++j) {
            int read_objects = fscanf(f, "%lf", &matrix->matr[i][j]);
            if (read_objects != 1) {
                free_matrix(matrix);
                fclose(f);
                fprintf(stderr, "Can't read element\n");
                return NULL;
            }
        }
    }
    fclose(f);
    return matrix;
}

Matrix* create_matrix(size_t rows, size_t cols) {
    Matrix* matrix = (Matrix*)malloc(sizeof(Matrix));
    if (!matrix) {
        return NULL;
    }

    matrix->rows = rows;
    matrix->cols = cols;
    matrix->matr = (double**)malloc(matrix->rows * sizeof(double*));
    if (!matrix->matr) {
        free(matrix);
        return NULL;
    }

    for (int i = 0; i < matrix->rows; ++i) {
        matrix->matr[i] = (double*)malloc(matrix->cols * sizeof(double));
        if (!matrix->matr[i]) {
            for (int k = 0; k < i; ++k) {
                free(matrix->matr[k]);
            }
            free(matrix->matr);
            free(matrix);
            return NULL;
        }
    }
    return matrix;
}

int free_matrix(Matrix* matrix) {
    for (int i = 0; i < matrix->rows; i++) {
        free(matrix->matr[i]);
    }
    free(matrix->matr);
    free(matrix);
    return 0;
}

int get_rows(const Matrix* matrix, size_t* rows) {
    if (!rows) {
        return 42;
    }
    *rows = matrix->rows;
    return 0;
}

int get_cols(const Matrix* matrix, size_t* cols) {
    if (!cols) {
        return 42;
    }
    *cols = matrix->cols;
    return 0;
}

int get_elem(const Matrix* matrix, int row, int col, double* val) {
    if (!val || row > matrix->rows || col > matrix->cols || row <  0 || col < 0) {
        return 42;
    }
    *val = matrix->matr[row][col];
    return 0;
}

int set_elem(Matrix* matrix, int row, int col, double val) {
    if (row > matrix->rows || col > matrix->cols || row <  0 || col < 0) {
        return 42;
    }
    matrix->matr[row][col] = val;
    return 0;
}

Matrix* mul_scalar(const Matrix* matrix, double val) {
    Matrix* new_matrix = create_matrix(matrix->rows, matrix->cols);
    if (!new_matrix) {
        return NULL;
    }

    for (int i = 0; i < new_matrix->rows; ++i) {
        for (int j = 0; j < new_matrix->cols; ++j) {
            new_matrix->matr[i][j] = matrix->matr[i][j] * val;
        }
    }
    return new_matrix;
}

Matrix* transp(const Matrix* matrix) {
    Matrix* new_matrix = create_matrix(matrix->cols, matrix->rows);
    if (!new_matrix) {
        return NULL;
    }

    for (int i = 0; i < new_matrix->rows; ++i) {
        for (int j = 0; j < new_matrix->cols; ++j) {
            new_matrix->matr[i][j] = matrix->matr[j][i];
        }
    }
    return new_matrix;
}

Matrix* sum(const Matrix* l, const Matrix* r) {
    if (l->rows != r->rows || l->cols != r->cols) {
        return NULL;
    }

    Matrix* new_matrix = create_matrix(l->rows, l->cols);
    if (!new_matrix) {
        return NULL;
    }

    for (int i = 0; i < new_matrix->rows; ++i) {
        for (int j = 0; j < new_matrix->cols; ++j) {
            new_matrix->matr[i][j] = l->matr[i][j] + r->matr[i][j];
        }
    }
    return new_matrix;
}

Matrix* sub(const Matrix* l, const Matrix* r) {
    if (l->rows != r->rows || l->cols != r->cols) {
        return NULL;
    }

    Matrix* new_matrix = create_matrix(l->rows, l->cols);
    if (!new_matrix) {
        return NULL;
    }

    for (int i = 0; i < new_matrix->rows; ++i) {
        for (int j = 0; j < new_matrix->cols; ++j) {
            new_matrix->matr[i][j] = l->matr[i][j] - r->matr[i][j];
        }
    }
    return new_matrix;
}

Matrix* mul(const Matrix* l, const Matrix* r) {
    if (l->cols != r->rows) {
        return NULL;
    }

    Matrix* new_matrix = create_matrix(l->rows, r->cols);
    if (!new_matrix) {
        return NULL;
    }

    for (int i = 0; i < new_matrix->rows; ++i) {
        for (int j = 0; j < new_matrix->cols; ++j) {
            double sum = 0;
            for (int k = 0; k < l->cols; ++k) {
                sum += l->matr[i][k] * r->matr[k][j];
            }
            new_matrix->matr[i][j] = sum;
        }
    }
    return new_matrix;
}

Matrix* get_minor(const Matrix* matrix, int row, int col) {
    Matrix* minor = create_matrix(matrix->rows - 1, matrix->cols - 1);
    if (!minor) {
        return NULL;
    }

    short row_flag = 0;
    for (int i = 0; i < minor->rows; ++i) {
        if (i == row) {
            row_flag = 1;
        }
        short col_flag = 0;
        for (int j = 0; j < minor->cols; ++j) {
            if (j == col) {
                col_flag = 1;
            }
            minor->matr[i][j] = matrix->matr[i+row_flag][j+col_flag];
        }
    }
    return minor;
}

int det(const Matrix* matrix, double* val) {
    if (matrix->rows != matrix->cols || !val) {
       return 42;
    }

    if (matrix->rows == 1) {
        * val = matrix->matr[0][0];
    } else if (matrix->rows == 2) {
        * val = matrix->matr[0][0] * matrix->matr[1][1] - (matrix->matr[1][0] * matrix->matr[0][1]);
    } else {
        double sum = 0.0;
        int koef = 1;

        for (int i = 0; i < matrix->rows; ++i) {
            Matrix *minor = get_minor(matrix, 0, i);
            if (!minor || det(minor, val)) {
                return 42;
            }

            sum += koef * (* val) * matrix->matr[0][i];
            free_matrix(minor);
            koef *= -1;
        }
        * val = sum;
    }
    return 0;
}

Matrix* adj(const Matrix* matrix) {
    if (matrix->rows != matrix->cols) {
       return NULL;
    }

    Matrix* new_matrix = create_matrix(matrix->rows, matrix->cols);
    if (!new_matrix) {
        return NULL;
    }

    for (int i = 0; i < new_matrix->rows; ++i) {
        for (int j = 0; j < new_matrix->cols; ++j) {
            Matrix *minor = get_minor(matrix, i, j);
            if (!minor) {
                free_matrix(new_matrix);
                return NULL;
            }

            double val = 0;
            if (det(minor, &val)) {
                free_matrix(minor);
                free_matrix(new_matrix);
                return NULL;
            }

            free_matrix(minor);
            new_matrix->matr[j][i] = pow(-1, i+j) * val;
        }
    }
    return new_matrix;
}

Matrix* inv(const Matrix* matrix) {
    if (matrix->rows != matrix->cols) {
       return NULL;
    }

    Matrix* new_matrix = create_matrix(matrix->rows, matrix->cols);
    if (!new_matrix) {
        return NULL;
    }

    if (new_matrix->rows == 1) {
        new_matrix->matr[0][0] = pow(matrix->matr[0][0], -1);
        return new_matrix;
    }

    double val = 0;
    det(matrix, &val);
    if (!val) {
        free_matrix(new_matrix);
        return NULL;
    }

    Matrix* adjugate = adj(matrix);
    if (!adjugate) {
        free_matrix(new_matrix);
        return NULL;
    }

    for (int i = 0; i < new_matrix->rows; ++i) {
        for (int j = 0; j < new_matrix->cols; ++j) {
            new_matrix->matr[i][j] = adjugate->matr[i][j] / val;
        }
    }

    free_matrix(adjugate);
    return new_matrix;
}
