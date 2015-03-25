%module BogusInterface

%include start.i
%include solverOptions.i
%import FCLib.i

%{
#include <fclib.h>
#include "bogus_interface.hpp"
%}
%include bogus_interface.hpp

