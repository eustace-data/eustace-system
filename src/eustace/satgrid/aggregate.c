/* aggregate.c: Python extension module for efficient aggregation */

/* required to avoid errors when including NumPy */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdio.h>

/* function declarations and corresponding documentation strings */

PyObject* aggregate_sum(PyObject* self, PyObject* args, PyObject* kw);
static char* aggregate_sum_keywords[] = { "bins", "values", "bin_indices", NULL };
static char aggregate_sum_docs[] = "aggregate_sum(bins, values, bin_indices): sum values into bins at specified bin indices.\n";

PyObject* aggregate_count_sum_min_max(PyObject* self, PyObject* args, PyObject* kw);
static char* aggregate_count_sum_min_max_keywords[] = { "bins_count", "bins_sum", "bins_min", "bins_max", "values", "bin_indices", NULL };
static char aggregate_count_sum_min_max_docs[] = "aggregate_count_sum_min_max(bins_count, bins_sum, bins_min, bins_max, values, bin_indices): put statistics of values into bins at specified bin indices.\n";

PyObject* aggregate_constant(PyObject* self, PyObject* args, PyObject* kw);
static char* aggregate_constant_keywords[] = { "bins", "k", "bin_indices", NULL };
static char aggregate_constant_docs[] = "aggregate_constant(bins, k, bin_indices): add constant k to bins at specified bin indices.\n";

PyObject* aggregate_count(PyObject* self, PyObject* args, PyObject* kw);
static char* aggregate_count_keywords[] = { "bins", "bin_indices", NULL };
static char aggregate_count_docs[] = "aggregate_count(bins, bin_indices): compute occupancy of specified bins.\n";

PyObject* aggregate_sum_dev_sq(PyObject* self, PyObject* args, PyObject* kw);
static char* aggregate_sum_dev_sq_keywords[] = { "bins", "bin_mean", "values", "bin_indices", NULL };
static char aggregate_sum_dev_sq_docs[] = "aggregate_sum_dev_sq(bins, bin_mean, values, bin_indices): compute sum square deviations about a known mean.\n";

/* index of all functions in this module*/
static PyMethodDef aggregate_funcs[] = {
  { "aggregate_sum", (PyCFunction)aggregate_sum, METH_VARARGS|METH_KEYWORDS, aggregate_sum_docs},
  { "aggregate_count_sum_min_max", (PyCFunction)aggregate_count_sum_min_max, METH_VARARGS|METH_KEYWORDS, aggregate_count_sum_min_max_docs},
  { "aggregate_constant", (PyCFunction)aggregate_constant, METH_VARARGS|METH_KEYWORDS, aggregate_constant_docs},
  { "aggregate_count", (PyCFunction)aggregate_count, METH_VARARGS|METH_KEYWORDS, aggregate_count_docs},
  { "aggregate_sum_dev_sq", (PyCFunction)aggregate_sum_dev_sq, METH_VARARGS|METH_KEYWORDS, aggregate_sum_dev_sq_docs},
  { NULL }
};

/* module initialisation */
void initaggregate(void)
{
  Py_InitModule3("aggregate", aggregate_funcs, "extension module for efficient aggregation");
  import_array();
}

PyObject* aggregate_sum(PyObject* self, PyObject* args, PyObject* kw)
{
  PyArrayObject* py_bins = NULL;
  PyArrayObject* py_values = NULL;
  PyArrayObject* py_bin_indices = NULL;
  npy_int32 num_indices = 0;
  npy_int32 loopcount = 0;
  npy_float32* bins = NULL;
  npy_float32* iter_value = NULL;
  npy_int32* iter_bin_index = NULL;
  
  if (!PyArg_ParseTupleAndKeywords(args, kw, "O!O!O!", aggregate_sum_keywords,
				   &PyArray_Type, &py_bins, 
				   &PyArray_Type, &py_values, 
				   &PyArray_Type, &py_bin_indices))
  {
    return NULL;
  }

  num_indices = PyArray_SIZE(py_bin_indices);
  bins = (npy_float32*)PyArray_DATA(py_bins);
  iter_value = (npy_float32*)PyArray_DATA(py_values);
  iter_bin_index = (npy_int32*)PyArray_DATA(py_bin_indices);

  for (loopcount = 0; loopcount < num_indices; loopcount++, iter_value++, iter_bin_index++)
    bins[*iter_bin_index] += *iter_value;

  Py_RETURN_NONE;
}

PyObject* aggregate_count_sum_min_max(PyObject* self, PyObject* args, PyObject* kw)
{
  PyArrayObject* py_bins_count = NULL;
  PyArrayObject* py_bins_sum = NULL;
  PyArrayObject* py_bins_min = NULL;
  PyArrayObject* py_bins_max = NULL;
  PyArrayObject* py_values = NULL;
  PyArrayObject* py_bin_indices = NULL;
  npy_int32 num_indices = 0;
  npy_int32 loopcount = 0;
  npy_int32* bins_count = NULL;
  npy_float32* bins_sum = NULL;
  npy_float32* bins_min = NULL;
  npy_float32* bins_max = NULL;
  npy_float32* iter_value = NULL;
  npy_int32* iter_bin_index = NULL;
  npy_float32* ptr_min = NULL;
  npy_float32* ptr_max = NULL;
  
  if (!PyArg_ParseTupleAndKeywords(args, kw, "O!O!O!O!O!O!", aggregate_count_sum_min_max_keywords,
				   &PyArray_Type, &py_bins_count, 
				   &PyArray_Type, &py_bins_sum, 
				   &PyArray_Type, &py_bins_min, 
				   &PyArray_Type, &py_bins_max, 
				   &PyArray_Type, &py_values, 
				   &PyArray_Type, &py_bin_indices))
  {
    return NULL;
  }

  num_indices = PyArray_SIZE(py_bin_indices);
  bins_count = (npy_int32*)PyArray_DATA(py_bins_count);
  bins_sum = (npy_float32*)PyArray_DATA(py_bins_sum);
  bins_min = (npy_float32*)PyArray_DATA(py_bins_min);
  bins_max = (npy_float32*)PyArray_DATA(py_bins_max);
  iter_value = (npy_float32*)PyArray_DATA(py_values);
  iter_bin_index = (npy_int32*)PyArray_DATA(py_bin_indices);

  for (loopcount = 0; loopcount < num_indices; loopcount++, iter_value++, iter_bin_index++)
  {
    bins_count[*iter_bin_index]++;
    bins_sum[*iter_bin_index] += *iter_value;
    ptr_min = bins_min + *iter_bin_index;
    ptr_max = bins_max + *iter_bin_index;
    if (*iter_value < *ptr_min)
      *ptr_min = *iter_value;
    if (*iter_value > *ptr_max)
      *ptr_max = *iter_value;
  }

  Py_RETURN_NONE;
}

PyObject* aggregate_constant(PyObject* self, PyObject* args, PyObject* kw)
{
  PyArrayObject* py_bins = NULL;
  PyArrayObject* py_bin_indices = NULL;
  npy_int32 num_indices = 0;
  npy_int32 loopcount = 0;
  npy_float32* bins = NULL;
  npy_float32 k = 0.0f;
  npy_int32* iter_bin_index = NULL;
  
  if (!PyArg_ParseTupleAndKeywords(args, kw, "O!fO!", aggregate_constant_keywords,
				   &PyArray_Type, &py_bins, 
				   &k, 
				   &PyArray_Type, &py_bin_indices))
  {
    return NULL;
  }

  num_indices = PyArray_SIZE(py_bin_indices);
  bins = (npy_float32*)PyArray_DATA(py_bins);
  iter_bin_index = (npy_int32*)PyArray_DATA(py_bin_indices);

  for (loopcount = 0; loopcount < num_indices; loopcount++, iter_bin_index++)
    bins[*iter_bin_index] += k;

  Py_RETURN_NONE;
}

PyObject* aggregate_count(PyObject* self, PyObject* args, PyObject* kw)
{
  PyArrayObject* py_bins = NULL;
  PyArrayObject* py_bin_indices = NULL;
  npy_int32 num_indices = 0;
  npy_int32 loopcount = 0;
  npy_int32* bins = NULL;
  npy_int32* iter_bin_index = NULL;
  
  if (!PyArg_ParseTupleAndKeywords(args, kw, "O!O!", aggregate_count_keywords,
				   &PyArray_Type, &py_bins, 
				   &PyArray_Type, &py_bin_indices))
  {
    return NULL;
  }

  num_indices = PyArray_SIZE(py_bin_indices);
  bins = (npy_int32*)PyArray_DATA(py_bins);
  iter_bin_index = (npy_int32*)PyArray_DATA(py_bin_indices);

  for (loopcount = 0; loopcount < num_indices; loopcount++, iter_bin_index++)
    bins[*iter_bin_index]++;

  Py_RETURN_NONE;
}

PyObject* aggregate_sum_dev_sq(PyObject* self, PyObject* args, PyObject* kw)
{
  PyArrayObject* py_bins = NULL;
  PyArrayObject* py_bin_mean = NULL;
  PyArrayObject* py_values = NULL;
  PyArrayObject* py_bin_indices = NULL;
  npy_int32 num_indices = 0;
  npy_int32 loopcount = 0;
  npy_float32* bins = NULL;
  npy_float32* bin_mean = NULL;
  npy_float32* iter_value = NULL;
  npy_int32* iter_bin_index = NULL;
  npy_float32 dev = 0.0f;
  
  if (!PyArg_ParseTupleAndKeywords(args, kw, "O!O!O!O!", aggregate_sum_dev_sq_keywords,
				   &PyArray_Type, &py_bins, 
				   &PyArray_Type, &py_bin_mean, 
				   &PyArray_Type, &py_values, 
				   &PyArray_Type, &py_bin_indices))
  {
    return NULL;
  }

  num_indices = PyArray_SIZE(py_bin_indices);
  bins = (npy_float32*)PyArray_DATA(py_bins);
  bin_mean = (npy_float32*)PyArray_DATA(py_bin_mean);
  iter_value = (npy_float32*)PyArray_DATA(py_values);
  iter_bin_index = (npy_int32*)PyArray_DATA(py_bin_indices);

  for (loopcount = 0; loopcount < num_indices; loopcount++, iter_value++, iter_bin_index++)
  {
    dev = *iter_value - bin_mean[*iter_bin_index];
    bins[*iter_bin_index] += dev * dev;
  }

  Py_RETURN_NONE;
}
