#!/usr/bin/env python3

import sys
import os
import subprocess
import pprint

from plot import plot

# --------------------------------- ADD YOUR BENCHMARKS HERE
BENCHMARKS = [{"name": "benchmark_test_macro", "configs":[("",""),("-r",{"1","2"})]},
         {"name": "benchmark_test_xyz_f", "configs":[("",""),("-r",{"2","3"})]},
         {"name": "benchmark_polyvest", "configs":[("",{""})]}
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
   results = []
   bname = benchmark["name"]
   bconfigs = benchmark["configs"]
   for config in bconfigs:
      for option in config[1]:
         print("# Running Benchmark '{}' with config '{} {}'...".format(bname, config[0], option));
         myenv = os.environ;
         #myenv["OMP_NUM_THREADS"] = str(nproc); # change env
         proc = subprocess.Popen([sys.path[0]+"/"+bname, config[0]+" "+option], stdout=subprocess.PIPE, stderr=subprocess.PIPE, env = myenv);
         f = open(sys.path[0]+"/out/"+bname+config[0]+"_"+option+".out", "w")
         for line in proc.stdout:
            try:
               dict = eval(line)
               results.append((("_"+config[0]+"_"+option).replace("__",""),dict))
               f.write(str(dict)+'\n')
            except:
               f.write(line.decode('utf-8'))
         f.close()
   return results

do_plot = True
plot_name = None
results = []

if len(sys.argv) == 2:
   plot_name = (DO_BENCHMARKS[0])["name"]

for benchmark in DO_BENCHMARKS:
   result = run_benchmark(benchmark)
   results.extend(result)

if do_plot and results:
   plot(sys.path[0], plot_name, results)