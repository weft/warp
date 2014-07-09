//wgeometry.i
%module wgeometry
%{ 
#include <vector>
#include <iostream>
#include <sstream>
#include <cmath>
#include <assert.h>
#include <time.h>
#include <string.h>
#include "datadef.h"
#include "primitive.h"
#include "wgeometry.h"
%}
%include "wgeometry.h"
%include "std_vector.i"
namespace std {
%template(Unsigned)   vector < unsigned >;
%template(Float)      vector < float >;
}
