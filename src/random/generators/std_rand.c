#include <stdlib.h>
#include <time.h>
#include "std_rand.h"


void std_init(void *seed){
    srand((unsigned) time(seed));
}

inline uint32_t std_rand(){
    return rand();
}
