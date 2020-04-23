#!/usr/bin/env python3

import sys
import os
import subprocess
import pprint
import operator

from plot import plot

# --------------------------------- ADD YOUR BENCHMARKS HERE

# --- Available functions for CLI
# List all function version choices for a benchmark that you want to execute here.
xyz_f = ["xyz_f1", "xyz_f2"]
for index, item in enumerate(xyz_f):
   xyz_f[index] = "xyz_f="+item

dotProduct = ["ref","auto1","auto2","vec1","2acc"]
for index, item in enumerate(dotProduct):
   dotProduct[index] = "dotProduct="+item


# --- Benchmarks
'''
	OUTDATED!!!
    id:            unique benchmark id
    name:          name of the file to benchmark
    run_configs:   CLI options that determine different functions to be executed
                   format: ("-option_name", {"option_value0", "option_value1", ...})
    input_configs: CLI options that determine different inputs for functions
                   format: ("-option_name", {{option_set0}, "option_value1", ...})
                           option_set can either be different option_values or ONE entry of the predefined structure {"i", operator, increment, lower_bound, upper_bound}.
                           In this case, there is a loop that iterates over the desired range and only 
'''
BENCHMARKS = [{"id": 0, "name": "benchmark_test_macro", "fun_configs":[], "run_configs":["r=10"], "input_configs":[]},
         {"id": 1, "name": "benchmark_test_xyz_f", "fun_configs": xyz_f, "run_configs":[], "input_configs":[]},
         {"id": 2, "name": "benchmark_polyvest", "fun_configs":[], "run_configs":[], "input_configs":[]},
         {"id": 3, "name": "benchmark_dotProduct", "fun_configs":dotProduct, "run_configs":["r=100000"], "input_configs":[("n",["i", operator.add, 32, 1, 50])]}
        ];

# --- Functions that should be compared

COMPARE = [2]


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
   compare_results = []
   bid = benchmark["id"]
   bname = benchmark["name"]
   bfconfigs = benchmark["fun_configs"]
   brconfigs = benchmark["run_configs"]
   biconfigs = benchmark["input_configs"]
   if not bfconfigs:
      bfconfigs.append("")
   if not brconfigs:
      brconfigs.append("")
   if not biconfigs:
      biconfigs.append(("",[""]))
   for fconfig in bfconfigs:
      for rconfig in brconfigs:
         for iconfigs in biconfigs:
            icname = iconfigs[0]
            ioptions = iconfigs[1]
            if ioptions[0] == "i":
               i = ioptions[3]
               while i <= ioptions[4]:
                  ival = str(i)
                  if fconfig == "":
                     config_string = "-b {}, {}={},".format(rconfig, icname, ival).replace(", =,", ",").replace(", -b ,", ",").replace("-b ,", "")
                  else:
                     config_string = "-f {}, -b {}, {}={},".format(fconfig, rconfig, icname, ival).replace(", =,", ",").replace(", -b ,", ",").replace("-b ,", "")
                  print("# Running Benchmark '{}' with config '{}'...".format(bname, config_string));
                  myenv = os.environ;
                  proc = subprocess.Popen([sys.path[0]+"/"+bname, config_string], stdout=subprocess.PIPE, stderr=subprocess.PIPE, env = myenv);
                  #proc = subprocess.Popen([sys.path[0]+"/"+bname, "-r=10"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, env = myenv);
                  f = open((sys.path[0]+"/out/"+bname+config_string.replace(" ", "_")).replace("__","")+".out", "w")
                  for line in proc.stdout:
                     try:
                        dict = eval(line)
                        if bid in COMPARE:
                           compare_results.append((("_"+config_string.replace(" ", "_")).replace("__",""),dict))
                        f.write(str(dict)+'\n')
                     except:
                        f.write(line.decode('utf-8'))
                  f.close()
                  i = ioptions[1](i, ioptions[2])
            else:
               for ival in ioptions:
                  if fconfig == "":
                     config_string = "-b {}, {}={},".format(rconfig, icname, ival).replace(", =,", ",").replace(", -b ,", ",").replace("-b ,", "")
                  else:
                     config_string = "-f {}, -b {}, {}={},".format(fconfig, rconfig, icname, ival).replace(", =,", ",").replace(", -b ,", ",").replace("-b ,", "")
                  print("# Running Benchmark '{}' with config '{}'...".format(bname, config_string));
                  myenv = os.environ;
                  proc = subprocess.Popen([sys.path[0]+"/"+bname, config_string], stdout=subprocess.PIPE, stderr=subprocess.PIPE, env = myenv);
                  f = open((sys.path[0]+"/out/"+bname+config_string.replace(" ", "_").replace(",", "")).replace("__","")+".out", "w")
                  for line in proc.stdout:
                     try:
                        dict = eval(line)
                        if bid in COMPARE:
                           compare_results.append(((config_string.replace(" ", "_").replace(",", "")).replace("__",""),dict))
                        f.write(str(dict)+'\n')
                     except:
                        f.write(line.decode('utf-8'))
                  f.close()
   return compare_results

do_plot = True
plot_name = None
compare_results = []

if len(sys.argv) == 2:
   plot_name = (DO_BENCHMARKS[0])["name"]

for benchmark in DO_BENCHMARKS:
   result = run_benchmark(benchmark)
   compare_results.extend(result)

if do_plot and compare_results:
   plot(sys.path[0], plot_name, compare_results)
