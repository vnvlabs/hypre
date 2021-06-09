#include <stddef.h>
#include "VnV.h"

#ifndef HYPRE_SEQUENTIAL
INJECTION_EXECUTABLE(VnVHypre, VNV, mpi) 
#else
INJECTION_EXECUTABLE(VnVHypre, VNV, serial) 
#endif

static const char* hypre_vnv_schema = "{\"type\": \"object\", \"required\":[]}";

/**
 * VnV allows users to set options in hypre using the input file. This callback
 * registers the json schema for the available options with the toolkit and defines
 * the options callback. Because we are in C, we are using the C interface. 
 * 
 * Life would be way easier if we can just compile this file with c++ :)
 *
 * TODO: Add options to the schema and parse them in this function.
 */ 
INJECTION_OPTIONS(VnVHypre, hypre_vnv_schema) {

}


/** For testing purposes only -- if you are seeing this, i forgot to remove it
 *
 */

int vnv_test_function(int x) {
   

  INJECTION_LOOP_BEGIN("VnVHypre", VWORLD("VnVHypre"), "SanityCheck", x)
  for (int i = 0; i < 10; i++) {
    x += i;
    INJECTION_LOOP_ITER("VnVHypre","SanityCheck", "inner");
  }

  INJECTION_LOOP_END("VnVHypre","SanityCheck");
  return x;
}
