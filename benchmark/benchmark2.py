#!/usr/bin/env python3

import sys
import os
import subprocess
import pprint
import operator
import itertools
import re

from plot import plot

# --------------------------------- ADD YOUR BENCHMARKS HERE

STD_REPS = "10"

# --- Available functions for CLI
# List all function version choices for a benchmark that you want to execute here.
xyz_f = ["xyz_f1", "xyz_f2"]
for index, item in enumerate(xyz_f):
   xyz_f[index] = "xyz_f="+item

dotProduct = ["ref","auto1","auto2","vec1","2acc"]
for index, item in enumerate(dotProduct):
   dotProduct[index] = "dotProduct="+item


intersectbodies = [ "cube_rot_r1.0_3", "cube_rot_r1.0_10", "cube_rot_r1.0_20",  "cube_rot_r1.0_40"]
intersectdims = {"cube_rot_r1.0_10": '10', "cube_rot_r1.0_3": '3', "cube_rot_r1.0_20": '20', "cube_rot_r1.0_40": '40'}

intersectEbodies = [ "ball_r1.0_3", "ball_r1.0_10", "ball_r1.0_20",  "ball_r1.0_40"]
intersectEdims = {"ball_r1.0_10": '10', "ball_r1.0_3": '3', "ball_r1.0_20": '20', "ball_r1.0_40": '40'}

# --- Benchmarks
'''
    id:            unique benchmark id, used for COMPARE
    name:          name of the file to benchmark
    config:        list of benchmark configs containing:

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
         xoption:        Labels the x-axis for benchmarks that have an input config
                        format: same as one element of the input config, i.e.
                                ("option-char", ["option-val0-printable, ..., option-valn-printable"])
                        note that length of list must match the one of the corresponding input config
                        we match by index!
'''
BENCHMARKS = [
   {"id": 0,
    "name": "benchmark_dotProduct",
    "config": [
       {
          "fun_configs": dotProduct,
          "run_configs": ["r=100000"],
          "input_configs": [("n", [2**i for i in range(0,7)])],
       }
    ],
    "xoption": ("n", {str(2**i): str(2**i) for i in range(0,7)}),
    "title": ["Runtime Comparison", "Performance comparison"],
    "xlabel": ["n", "n"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)"]
   },
   {"id": 1,
    "name": "benchmark_intersect",
    "config": [       
       {
          "fun_configs": ["Polytope_intersectCoord=cached_ref", "Polytope_intersectCoord=ref"],
          "run_configs": ["intersect=intersectCoord,polytopeTranspose=false"],
          "input_configs": [("generator", intersectbodies)]
       },
       {
          "fun_configs": ["PolytopeT_intersectCoord=cached_nc1", "PolytopeT_intersectCoord=ref"],
          "run_configs": ["intersect=intersectCoord,polytopeTranspose=true"],
          "input_configs": [("generator", intersectbodies)]
       },
       {
          "fun_configs": [],
          "run_configs": ["intersect=intersect,polytopeTranspose=false", "intersect=intersect,polytopeTranspose=true"],
          "input_configs": [("generator", intersectbodies)]
       }
    ],
    "xoption": ("generator", intersectdims),
    "title": ["Runtime Comparison", "Performance comparison"],
    "xlabel": ["dim", "dim"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)"]
    #"xoption": ("generator", intersectEdims)
   }
]


assert(len(set(map(lambda t: t["name"], BENCHMARKS))) == len(BENCHMARKS) and "benchmarks don't have unique names!")


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


         
"""
config contains the following:
fun_config :: [string]
run_config :: [string]
input_conifg :: [(char, [(readable, char)])]

from this we create the following list of config strings
'-f "<fun>" -n "<run,c=str(cconf)>"' 
for 
   fun : fun_config, 
   run : run_config,
   (c, conf) : input_config
   cconf : conf
"""
def get_config(config):
   funs = [('-f', fun) for fun in config["fun_configs"]] or [('','')]
   inputs = [('-b', ','.join(c_prod)) for c_prod in
             itertools.product(
                *(list(map(
                      lambda c: ['{}={}'.format(c[0], conf) for conf in c[1]],
                      config["input_configs"]
                )) + [config["run_configs"]])
             )
   ] or [('','')]
   return itertools.product(funs, inputs)


"""
benchmark contains 
config :: [dict]

we extract all config strings from the elements of this list and concat them
"""
def get_configs(benchmark):
   return [conf_string
           for conf in benchmark["config"]
           for conf_string in get_config(conf)
   ]


"""
get x-label from benchmark_string by:
- matching on the option tag
- getting the value
- looking up the label of the value in the xoption

this is ugly, but should work if we only choose labels on input_config...  
"""
def get_label(xoption, benchmark_string):
   pattern = '\-b\s*".*{}=([^,"]*)'.format(xoption[0])
   strval = re.search(pattern, benchmark_string).group(1)
   return xoption[1][strval]
   


"""
expects a list of strings as created by get_configs
this allows concatenating benchmarks -> get_configs creates product of config options, concatenating benchmark strings creates sum of config options
"""
def run_benchmark(bid, bname, config_strings):
   results = []

   for runconf,inputconf in config_strings:
      config_string = '{} "{}" {} "{}"'.format(runconf[0], runconf[1], inputconf[0], inputconf[1])
      config_string_printable = (
         config_string
         .replace(" ", "_")
         .replace(",", "")
      ).replace("__","")

      print('# Running Benchmark \n{} {}\n'.format(sys.path[0]+"/"+bname, config_string));
      myenv = os.environ;
      proc = subprocess.Popen(
         [sys.path[0]+"/"+bname,
          runconf[0],
          runconf[1],
          inputconf[0],
          inputconf[1]],
         stdout=subprocess.PIPE,
         #stderr=subprocess.PIPE,
         env = myenv
      );
      f = open(sys.path[0]+ "/out/" + bname + config_string_printable + ".out", "w")
      f.write(config_string + "\n")
      for line in proc.stdout:
         print(line)
         try:
            dict = eval(line)
            results.append((config_string, dict))
            f.write(str(dict)+'\n')
         except:
            f.write(line.decode('utf-8'))
      f.close()

   return results


for benchmark in DO_BENCHMARKS:
   bname = benchmark["name"]
   result = run_benchmark(benchmark["id"],
                          bname,
                          get_configs(benchmark)
   )
   # get x-axis labels and add them to data 
   result = list(map(lambda res: (*res, get_label(benchmark["xoption"], res[0])), result))
   plot(sys.path[0], bname, result, benchmark["xoption"][0], benchmark["title"], benchmark["xlabel"], benchmark["ylabel"])


'''
       ,
       {
          "fun_configs": ["Ellipsoid_intersectCoord=cached_ref", "Ellipsoid_intersectCoord=cached_reord"],
          "run_configs": ["r=100000,intersect=intersectCoord"],
          "input_configs": [("generator", intersectEbodies)]
       }
'''
