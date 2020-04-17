#!/usr/bin/env python3

import sys
import os
import subprocess
import pprint

# --------------------------------- ADD YOUR BENCHMARKS HERE
BENCHMARKS = [{"name": "benchmark_test_macro", "configs":["","-r 2"]},
         {"name": "benchmark_test_xyz_f", "configs":["","-r 2"]},
         {"name": "benchmark_polyvest", "configs":[""]}
        ];


# ---------------------------- Parse python args
DO_BENCHMARKS = BENCHMARKS;
if(len(sys.argv)>1):
   DO_BENCHMARKS = [];
   
   for i in range(1,len(sys.argv)):
      benchmark = sys.argv[i]
      upd = list(filter(lambda b: b["name"] == benchmark, BENCHMARKS))
      DO_BENCHMARKS += upd
      if len(upd) == 0:
         print("ERROR: benchmark '{}' is not available!".format(benchmark));
         print("  list of available benchmarks:");
         for b in BENCHMARKS:
            print("     ",b);
         sys.exit(1);

# ------ iterate over benchmarks chosen in DO_BENCHMARKS, subset of BENCHMARKS

def run_benchmark(benchmark):
   bname = benchmark["name"]
   bconfigs = benchmark["configs"]
   for config in bconfigs:
      print("# Running Benchmark '{}' with config '{}'...".format(bname, config));
      myenv = os.environ;
      #myenv["OMP_NUM_THREADS"] = str(nproc); # change env
      proc = subprocess.Popen([sys.path[0]+"/"+bname, config], stdout=subprocess.PIPE, stderr=subprocess.PIPE, env = myenv);
      f = open(sys.path[0]+"/out/"+bname+config.replace(' ', '_')+".out", "w")
      for line in proc.stdout:
         dict = eval(line)
         f.write(str(dict)+'\n')
      f.close()


for benchmark in DO_BENCHMARKS:
   run_benchmark(benchmark)

