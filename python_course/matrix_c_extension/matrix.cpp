#define PY_SSIZE_T_CLEAN
#include <Python.h>

extern PyTypeObject matrix_type;

typedef struct {
    PyObject_HEAD
    long int **matr;
    Py_ssize_t rows;
    Py_ssize_t cols;
} matrix_t;

static PyObject *
matrix_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    matrix_t *self;

    self = (matrix_t *)type->tp_alloc(type, 0);
    if (self == NULL) {
        return PyErr_NoMemory();
    }

    self->matr = NULL;
    self->rows = 0;
    self->cols = 0;
    return (PyObject *)self;
}

static int
matrix_init(matrix_t *self, PyObject *args, PyObject *kwds)
{
    PyObject *pList = NULL;

    if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &pList)) {
        PyErr_SetString(PyExc_TypeError, "argument must be a list of lists");
        return 1;
    }

    self->rows = PyList_Size(pList);
    if (self->rows < 1) {
        PyErr_SetString(PyExc_ValueError, "number of rows must be 1 at least");
        return 1;
    }

    self->matr = new long int *[self->rows];
    if (self->matr == NULL) {
        PyErr_SetString(PyExc_MemoryError, "no memory");
        return -1;
    }

    PyObject *pRow = NULL;
    PyObject *pItem = NULL;

    for (int i = 0; i < self->rows; ++i) {
        pRow = PyList_GetItem(pList, i);
        if(!PyList_Check(pRow)) {
            PyErr_SetString(PyExc_TypeError, "list items must be a list");
            return 1;
        }

        if (i == 0) {
            self->cols = PyList_Size(pRow);
            if (self->cols < 1) {
                PyErr_SetString(PyExc_ValueError, "number of cols must be 1 " \
                                                  "at least");
                return 1;
            }
        } else {
            Py_ssize_t n = PyList_Size(pRow);
            if (n != self->cols) {
                PyErr_SetString(PyExc_ValueError, "list items must be the " \
                                                  "same length");
                return 1;
            }
        }

        self->matr[i] = new long int[self->cols];
        if (self->matr[i] == NULL) {
            PyErr_SetString(PyExc_MemoryError, "no memory");
            return -1;
        }

        for (int j = 0; j < self->cols; ++j) {
            pItem = PyList_GetItem(pRow, j);
            if(!PyLong_Check(pItem)) {
                PyErr_SetString(PyExc_TypeError, "row-list items must be " \
                                                 "integers");
                return 1;
            }

            self->matr[i][j] = PyLong_AsLong(pItem);
        }
    }

    return 0;
}

static void
matrix_dealloc(matrix_t *self)
{
    if (self->matr != NULL) {
        for (int i = 0; i < self->rows; ++i) {
            delete [] self->matr[i];
        }
        delete [] self->matr;

        self->matr = NULL;
    }

    Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject *
matrix_repr(matrix_t *self)
{
    _PyUnicodeWriter writer;
    _PyUnicodeWriter_Init(&writer);
    writer.overallocate = 1;
    writer.min_length = 12;
    PyObject *el_str = NULL;

    if (_PyUnicodeWriter_WriteASCIIString(&writer, "<matrix { ", 10) < 0) {
        goto error;
    }

    for (int i = 0; i < self->rows; ++i) {
        if (_PyUnicodeWriter_WriteASCIIString(&writer, "[ ", 2) < 0) {
            goto error;
        }

        for (int j = 0; j < self->cols; ++j) {
            el_str = PyUnicode_FromFormat("%d ", self->matr[i][j]);
            if (_PyUnicodeWriter_WriteStr(&writer, el_str)) {
                Py_DECREF(el_str);
                goto error;
            }

            Py_DECREF(el_str);
        }

        if (_PyUnicodeWriter_WriteASCIIString(&writer, "] ", 2) < 0) {
            goto error;
        }
    }

    writer.overallocate = 0;
    if (_PyUnicodeWriter_WriteASCIIString(&writer, "}>", 2) < 0) {
        goto error;
    }

    return _PyUnicodeWriter_Finish(&writer);

error:
    PyErr_SetString(PyExc_UnicodeError, "write failed");
    Py_XDECREF(el_str);
    _PyUnicodeWriter_Dealloc(&writer);
    return NULL;
}

static int
matrix_traverse(matrix_t *self, visitproc visit, void *arg)
{
    return 0;
}

static matrix_t *
create_matrix(Py_ssize_t rows, Py_ssize_t cols)
{
    matrix_t *matrix = (matrix_t *)matrix_new(&matrix_type, NULL, NULL);
    matrix->rows = rows;
    matrix->cols = cols;
    matrix->matr = new long int *[rows];
    if (matrix->matr == NULL) {
        PyErr_SetString(PyExc_MemoryError, "no memory");
        return NULL;
    }

    for (int i = 0; i < rows; ++i) {
        matrix->matr[i] = new long int[cols];
        if (matrix->matr[i] == NULL) {
            PyErr_SetString(PyExc_MemoryError, "no memory");
            return NULL;
        }
    }

    return matrix;
}

static PyObject *
matrix_sum(matrix_t *self, PyObject *args)
{
    matrix_t *other = NULL;

    if (!PyArg_ParseTuple(args, "O!", &matrix_type, &other)) {
        PyErr_SetString(PyExc_TypeError, "argument must be matrix");
        return NULL;
    }

    matrix_t *matrix = create_matrix(self->rows, self->cols);

    for (int i = 0; i < self->rows; ++i) {
        for (int j = 0; j < self->cols; ++j) {
            matrix->matr[i][j] = self->matr[i][j] + other->matr[i][j];
        }
    }

    return (PyObject *)matrix;
}

static PyObject *
matrix_mul_scalar(matrix_t *self, PyObject *args)
{
    long int num = 0;

    if (!PyArg_ParseTuple(args, "l", &num)) {
        PyErr_SetString(PyExc_TypeError, "argument must be integer");
        return NULL;
    }

    matrix_t *matrix = create_matrix(self->rows, self->cols);

    for (int i = 0; i < self->rows; ++i) {
        for (int j = 0; j < self->cols; ++j) {
            matrix->matr[i][j] = self->matr[i][j] * num;
        }
    }

    return (PyObject *)matrix;
}

static PyObject *
matrix_div(matrix_t *self, PyObject *args)
{
    long int num = 0;

    if (!PyArg_ParseTuple(args, "l", &num)) {
        PyErr_SetString(PyExc_TypeError, "argument must be integer");
        return NULL;
    }

    if (num == 0) {
        PyErr_SetString(PyExc_ZeroDivisionError, "argument should not be 0");
        return NULL;
    }

    matrix_t *matrix = create_matrix(self->rows, self->cols);

    for (int i = 0; i < self->rows; ++i) {
        for (int j = 0; j < self->cols; ++j) {
            matrix->matr[i][j] = self->matr[i][j] / num;
        }
    }

    return (PyObject *)matrix;
}

static PyObject *
matrix_mul(matrix_t *self, PyObject *args)
{
    matrix_t *other = NULL;

    if (!PyArg_ParseTuple(args, "O!", &matrix_type, &other)) {
        PyErr_SetString(PyExc_TypeError, "argument must be matrix");
        return NULL;
    }

    if (self->cols != other->rows) {
        PyErr_SetString(PyExc_ValueError, "col numer of left matrix mast be" \
                        "the same as row number of right one");
        return NULL;
    }

    matrix_t *matrix = create_matrix(self->rows, other->cols);

    for (int i = 0; i < matrix->rows; ++i) {
        for (int j = 0; j < matrix->cols; ++j) {
            long int sum = 0;
            for (int k = 0; k < self->cols; ++k) {
                sum += self->matr[i][k] * other->matr[k][j];
            }
            matrix->matr[i][j] = sum;
        }
    }

    return (PyObject *)matrix;
}

static PyObject *
matrix_transp(matrix_t *self)
{
    matrix_t *matrix = create_matrix(self->cols, self->rows);

    for (int i = 0; i < matrix->rows; ++i) {
        for (int j = 0; j < matrix->cols; ++j) {
            matrix->matr[i][j] = self->matr[j][i];
        }
    }

    return (PyObject *)matrix;
}

static PyObject *
matrix_get(matrix_t *self, PyObject *args)
{
    Py_ssize_t row = 0;
    Py_ssize_t col = 0;

    if (!PyArg_ParseTuple(args, "nn", &row, &col)) {
        PyErr_SetString(PyExc_TypeError, "argument must be like: (row,col)");
        return NULL;
    }

    if (row < 0 || row > self->rows-1) {
        PyErr_SetString(PyExc_IndexError, "incorrect row index");
        return NULL;
    }
    if (col < 0 || col > self->cols-1) {
        PyErr_SetString(PyExc_IndexError, "incorrect col index");
        return NULL;
    }

    return PyLong_FromLong(self->matr[row][col]);
}

static int
matrix_contains(matrix_t *self, PyObject *arg)
{
    if (!PyLong_Check(arg)) {
        PyErr_SetString(PyExc_TypeError, "argument must be integer");
        return -1;
    }

    long int num = PyLong_AsLong(arg);

    for (int i = 0; i < self->rows; ++i) {
        for (int j = 0; j < self->cols; ++j) {
            if (self->matr[i][j] == num) {
                return 1;
            }
        }
    }

    return 0;
}

static PySequenceMethods matrix_as_sequence = {
    0,                                               /* sq_length */
    0,                                               /* sq_concat */
    0,                                               /* sq_repeat */
    0,                                               /* sq_item */
    0,                                               /* sq_slice */
    0,                                               /* sq_ass_item */
    0,                                               /* sq_ass_slice */
    (objobjproc)matrix_contains,                     /* sq_contains */
};

static PyMappingMethods matrix_as_mapping = {
    0,                                               /* mp_length */
    (binaryfunc)matrix_get,                          /* mp_subscript */
    0,                                               /* mp_ass_subscript */
};

static PyMethodDef matrix_methods[] = {
    {"sum", (PyCFunction)matrix_sum, METH_VARARGS},
    {"mul_scalar", (PyCFunction)matrix_mul_scalar, METH_VARARGS},
    {"div", (PyCFunction)matrix_div, METH_VARARGS},
    {"mul", (PyCFunction)matrix_mul, METH_VARARGS},
    {"transp", (PyCFunction)matrix_transp, METH_NOARGS},
    {NULL, NULL, 0, NULL}
};

PyTypeObject matrix_type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "matrix.Matrix",                                 /* tp_name */
    sizeof(matrix_t),                                /* tp_basic_size */
    0,                                               /* tp_itemsize */
    (destructor)matrix_dealloc,                      /* tp_dealloc */
    0,                                               /* tp_print */
    0,                                               /* tp_getattr */
    0,                                               /* tp_setattr */
    0,                                               /* tp_reserved */
    (reprfunc)matrix_repr,                           /* tp_repr */
    0,                                               /* tp_as_number */
    &matrix_as_sequence,                             /* tp_as_sequence */
    &matrix_as_mapping,                              /* tp_as_mapping */
    0,                                               /* tp_hash */
    0,                                               /* tp_call */
    0,                                               /* tp_str */
    0,                                               /* tp_getattro */
    0,                                               /* tp_setattro */
    0,                                               /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /* tp_flags */
    0,                                               /* tp_doc */
    (traverseproc)matrix_traverse,                   /* tp_traverse */
    0,                                               /* tp_clear */
    0,                                               /* tp_richcompare */
    0,                                               /* tp_weaklistoffset */
    0,                                               /* tp_iter */
    0,                                               /* tp_iternext */
    matrix_methods,                                  /* tp_methods */
    0,                                               /* tp_members */
    0,                                               /* tp_getset */
    0,                                               /* tp_base */
    0,                                               /* tp_dict */
    0,                                               /* tp_descr_get */
    0,                                               /* tp_descr_set */
    0,                                               /* tp_dictoffset */
    (initproc)matrix_init,                           /* tp_init */
    0,                                               /* tp_alloc */
    matrix_new,                                      /* tp_new */
    0,                                               /* tp_free */
};

static struct PyModuleDef matrix_module = {
    PyModuleDef_HEAD_INIT,
    "matrix",
    NULL,
    -1,
};

PyMODINIT_FUNC
PyInit_matrix(void)
{
    PyObject *m;
    if (PyType_Ready(&matrix_type) < 0)
        return NULL;

    m = PyModule_Create(&matrix_module);
    if (m == NULL)
        return NULL;

    Py_INCREF(&matrix_type);
    PyModule_AddObject(m, "Matrix", (PyObject *) &matrix_type);
    return m;
}
