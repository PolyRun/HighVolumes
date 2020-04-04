#include <iostream>

extern "C" { // must be included C stlye
#include "../../src/poly/volume.h"
#include "../../src/poly/preprocess.h"
}

#include "../../src/poly/volume_helper.hpp"

#include "../../src/util/cli.hpp"
#include "../../src/util/cli_functions.hpp"

int main(int argc, char** argv) {
   CLI cli(argc,argv,"test_volume_basics");
   CLIFunctionsVolume cliFun(cli);
  
   cliFun.preParse();
   if (!cli.parse()) {return -1;}
   cliFun.postParse();
   
   // -------------------------------- start tests


   // -------------------------------- end tests

}

