#include <stdlib.h>
#include "std_rand.h"


void std_init(int seed){
    srand((unsigned) time(seed));
}

uint32_t std_random(){
    return rand();
}
