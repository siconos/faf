%module BogusInterface

%include start.i
%include solverOptions.i
%import FCLib.i

%typemap(in) (double *data)
  (PyArrayObject* array =NULL, int is_new_object=0)
{
  array = obj_to_array_allow_conversion($input, NPY_DOUBLE,
                                        &is_new_object);
  if (!array
      || !require_native(array))
    SWIG_fail;

  $1 = (double *) array_data(array);
}

%apply (double *data) { (double *reaction) };
%apply (double *data) { (double *velocity) };
%apply (double *data) { (double *reactions) };
%apply (double *data) { (double *velocities) };

// 1 : numinputs=0 mandatory to avoid arg
%typemap(in, numinputs=0) (double *out3) (double result[3])
{
}

// 2 : check must be done after in
%typemap(check) (double *out3)
{
  $1 = result$argnum;
}

// 3 : return arg
%typemap(argout) (double *out3)
{
  $result = SWIG_Python_AppendOutput($result,PyFloat_FromDouble(result$argnum[0]));
  $result = SWIG_Python_AppendOutput($result,PyFloat_FromDouble(result$argnum[1]));
  $result = SWIG_Python_AppendOutput($result,PyFloat_FromDouble(result$argnum[2]));
}

// 1 : numinputs=0 mandatory to avoid arg
%typemap(in, numinputs=0) (double *out9) (double result[9])
{
}

// 2 : check must be done after in
%typemap(check) (double *out9)
{
  $1 = result$argnum;
}

// 3 : return arg
%typemap(argout) (double *out9)
{
  $result = SWIG_Python_AppendOutput($result,PyFloat_FromDouble(result$argnum[0]));
  $result = SWIG_Python_AppendOutput($result,PyFloat_FromDouble(result$argnum[1]));
  $result = SWIG_Python_AppendOutput($result,PyFloat_FromDouble(result$argnum[2]));

  $result = SWIG_Python_AppendOutput($result,PyFloat_FromDouble(result$argnum[3]));
  $result = SWIG_Python_AppendOutput($result,PyFloat_FromDouble(result$argnum[4]));
  $result = SWIG_Python_AppendOutput($result,PyFloat_FromDouble(result$argnum[5]));

  $result = SWIG_Python_AppendOutput($result,PyFloat_FromDouble(result$argnum[6]));
  $result = SWIG_Python_AppendOutput($result,PyFloat_FromDouble(result$argnum[7]));
  $result = SWIG_Python_AppendOutput($result,PyFloat_FromDouble(result$argnum[8]));
}

%apply (double *out9) { (double *out9_1) };
%apply (double *out9) { (double *out9_2) };

%{
#include <fclib.h>
#include "bogus_interface.hpp"
%}
%include bogus_interface.hpp



