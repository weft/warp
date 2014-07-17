//warp.i
%module warp
%{ 
#include <vector>
#include <iostream>
#include <sstream>
#include <cmath>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <cudpp_hash.h>
#include <curand.h>
#include "datadef.h"
#include "primitive.h"
#include "wgeometry.h"
#include "whistory.h"
%}
%include "wgeometry.h"
%include "whistory.h"
%include "std_vector.i"
%include "std_string.i"
namespace std {
%template(Unsigned)   vector < unsigned >;
%template(Float)      vector < float >;  
}
