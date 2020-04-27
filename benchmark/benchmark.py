#!/usr/bin/env python3

import sys
import os
import subprocess
import pprint
import operator

from plot import plot, plot_input

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
    id:            unique benchmark id, used for COMPARE
    name:          name of the file to benchmark
    fun_configs:   CLI options that determine the function to be selected
                   format: funListName -> Definded above
    run_configs:   CLI options that determine different benchmark configurations to be executed
                   format: ["opt-char0=opt-value0", "opt-char1=opt-value1", ...]
    input_configs: CLI options that determine different inputs for functions
                   format: ("option-char", ["option-val0", "option-val1", ...])
				           or the following format:
                           ("option-char", ["i", operator, increment, lower_bound, upper_bound].
                           In this case, there is a loop that iterates over the desired range
						Only one option-char is available
'''
BENCHMARKS = [{"id": 0, "name": "benchmark_test_macro", "fun_configs":[], "run_configs":["r=10"], "input_configs":[]},
         {"id": 1, "name": "benchmark_test_xyz_f", "fun_configs": xyz_f, "run_configs":[], "input_configs":[]},
         {"id": 2, "name": "benchmark_polyvest", "fun_configs":[], "run_configs":[], "input_configs":[]},
         {"id": 3, "name": "benchmark_dotProduct", "fun_configs":dotProduct, "run_configs":["r=1000000"], "input_configs":[("n",["i", operator.add, 10, 1, 50])]}
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
   results = []
   no_biconfig = False
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
      no_biconfig = True
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
                  if len(config_string) > 0:
                     while config_string[-1] == ",":
                        config_string = config_string[:-1]
                  print("# Running Benchmark '{}' with config '{}'...".format(bname, config_string));
                  myenv = os.environ;
                  proc = subprocess.Popen([sys.path[0]+"/"+bname, config_string], stdout=subprocess.PIPE, stderr=subprocess.PIPE, env = myenv);
                  f = open((sys.path[0]+"/out/"+bname+config_string.replace(" ", "_")).replace("__","")+".out", "w")
                  for line in proc.stdout:
                     try:
                        dict = eval(line)
                        results.append((fconfig, dict, ival))
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
                  if len(config_string) > 0:
                     while config_string[-1] == ",":
                        config_string = config_string[:-1]
                  print("# Running Benchmark '{}' with config '{}'...".format(bname, config_string));
                  myenv = os.environ;
                  proc = subprocess.Popen([sys.path[0]+"/"+bname, config_string], stdout=subprocess.PIPE, stderr=subprocess.PIPE, env = myenv);
                  f = open((sys.path[0]+"/out/"+bname+config_string.replace(" ", "_").replace(",", "")).replace("__","")+".out", "w")
                  for line in proc.stdout:
                     try:
                        dict = eval(line)
                        results.append(((config_string.replace(" ", "_").replace(",", "")).replace("__",""),dict))
                        f.write(str(dict)+'\n')
                     except:
                        f.write(line.decode('utf-8'))
                  f.close()
   
   if no_biconfig:
      plot(sys.path[0], bname, results)
   else:
      plot_input(sys.path[0], bname, results, icname)
   
   if bid in COMPARE:
      return results
   else:
      return []

do_plot = True
plot_name = "benchmark_comparison"
compare_results = []

for benchmark in DO_BENCHMARKS:
   result = run_benchmark(benchmark)
   compare_results.extend(result)

if do_plot and compare_results:
   plot(sys.path[0], plot_name, compare_results)
