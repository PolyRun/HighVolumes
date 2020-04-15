#!/usr/bin/env python3

import sys
import os
import subprocess
import pprint

# --------------------------------- ADD YOUR BENCHMARKS HERE
BENCHMARKS = [{"name": "benchmark_test_macro"},
         {"name": "benchmark_test_xyz_f"},
         {"name": "benchmark_polyvest"}
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
   print("# Running Benchmark '{}'...".format(bname));
   myenv = os.environ;
   #myenv["OMP_NUM_THREADS"] = str(nproc); # change env
   proc = subprocess.Popen((sys.path[0]+"/"+bname,), stdout=subprocess.PIPE, stderr=subprocess.PIPE, env = myenv);
   f = open(sys.path[0]+"/out/"+bname+".out", "w")
   for line in proc.stdout:
      f.write(line.decode('utf-8'))
   f.close()


for benchmark in DO_BENCHMARKS:
   run_benchmark(benchmark)

