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

# --- Available functions for CLI
# List all function version choices for a benchmark that you want to execute here.

dotProduct = ["ref","auto1","auto2","vec1","2acc"]
for index, item in enumerate(dotProduct):
   dotProduct[index] = "dotProduct="+item


intersectbodies = [ "cube_rot_r1.0_3", "cube_rot_r1.0_10", "cube_rot_r1.0_20",  "cube_rot_r1.0_40", "cube_rot_r1.0_60", "cube_rot_r1.0_100"]
intersectdims = {"cube_rot_r1.0_10": '10', "cube_rot_r1.0_3": '3', "cube_rot_r1.0_20": '20', "cube_rot_r1.0_40": '40', "cube_rot_r1.0_60": '60', "cube_rot_r1.0_100": '100'}


mbintersectbodies = ["cube_rot_r1.0_10", "cube_rot_r1.0_20",  "cube_rot_r1.0_40", "cube_rot_r1.0_60", "cube_rot_r1.0_100"]
mbintersectdims = {"cube_rot_r1.0_10": '10', "cube_rot_r1.0_20": '20', "cube_rot_r1.0_40": '40', "cube_rot_r1.0_60": '60', "cube_rot_r1.0_100": '100'}

intersectEbodies = [ "ball_r1.0_3", "ball_r1.0_10", "ball_r1.0_20", "ball_r1.0_40", "ball_r1.0_60", "ball_r1.0_100", "ball_r1.0_150", "ball_r1.0_200"]
intersectEdims = {"ball_r1.0_10": '10', "ball_r1.0_3": '3', "ball_r1.0_20": '20', "ball_r1.0_40": '40', "ball_r1.0_60": '60', "ball_r1.0_100": '100', "ball_r1.0_150": '150', "ball_r1.0_200": '200'}


intersectSparseDims = [4,5,10,20,40,60,100,150,200]
intersectSparseDims = {"2var_"+str(i):str(i) for i in intersectSparseDims}
intersectSparseBodies = [ name for (name,i) in intersectSparseDims.items()]

intersectSparse4Dims = [4,5,10,20,40,60,100,150,200]
intersectSparse4Dims = {"4var_"+str(i):str(i) for i in intersectSparse4Dims}
intersectSparse4Bodies = [ name for (name,i) in intersectSparse4Dims.items()]

cubeRotDims = [3,10,20,40,60,80,100]
cubeRotDims = {"cube_rot_r1.0_"+str(i):str(i) for i in cubeRotDims}
cubeRotBodies = [ name for (name,i) in cubeRotDims.items()]

jitTest = [4*i for i in range(1,20)]
jitTestDims = {str(i):str(4*i) for i in jitTest}
jitTest = [str(i) for i in jitTest]

randValTypes = ["random_int", "random_int_in_range", "random_double", "random_double_0_1", "random_double_normal", "random_double_in_range"]
randValTypes_ = {"random_int": '0', "random_int_in_range": '1', "random_double": '2', "random_double_0_1": '3', "random_double_normal": '4', "random_double_in_range": '5'}

rd_0_1 = ["ref","fast"]
for index, item in enumerate(rd_0_1):
   rd_0_1[index] = "rd_0_1="+item

# --- Benchmarks
'''
    name:          name of the benchmark, has to be unique
    executable:    name of the file (executable) that is used for this benchmark
    config:        list of benchmark configs containing:

         const_configs: CLI options that determine the algorithm constants to be selected
                        format: ["opt-char0=opt-value0", "opt-char1=opt-value1", ...]
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
         xoption:       Labels the ticks of the x-axis for benchmarks that have an input config
                        format: same as one element of the input config, i.e.
                                ("option-char", ["option-val0-printable, ..., option-valn-printable"])
                        note that length of list must match the one of the corresponding input config
                        we match by index!
         title:         Titles of the plots
         xlabel:        Labels the x-axis of the plots
         ylabel:        Labels the y-axis of the plots
'''
BENCHMARKS = [
   {"name": "benchmark_dotProduct",
    "executable": "benchmark_dotProduct",
    "config": [
       {
          "const_configs": [],
          "fun_configs": dotProduct,
          "run_configs": ["r=100000"],
          "input_configs": [("n", [2**i for i in range(0,7)])],
       }
    ],
    "xoption": ("n", {str(2**i): str(2**i) for i in range(0,7)}),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison"],
    "xlabel": ["n", "n", "n"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)"]
   },

   {"name": "intersect_polytope",
    "executable": "benchmark_intersect",
    "config": [       
       {
          "const_configs": [],
          "fun_configs": ["Polytope_intersectCoord=cached_ref", "Polytope_intersectCoord=ref"],
          "run_configs": ["intersect=intersectCoord,polytopeTranspose=false"],
          "input_configs": [("generator", intersectbodies)]
       },
       {
          "const_configs": [],
          #"fun_configs": ["PolytopeT_intersectCoord=cached_nc1", "PolytopeT_intersectCoord=ref", "PolytopeT_intersectCoord=cached_b_vec"],
          "fun_configs": ["PolytopeT_intersectCoord=cached_nc1", "PolytopeT_intersectCoord=cached_b_vec"],
          "run_configs": ["intersect=intersectCoord,polytopeTranspose=true"],
          "input_configs": [("generator", intersectbodies)]
       },
       #{
       #   "const_configs": [],
       #   "fun_configs": [],
       #   "run_configs": ["intersect=intersect,polytopeTranspose=false", "intersect=intersect,polytopeTranspose=true"],
       #   "input_configs": [("generator", intersectbodies)]
       #}
    ],
    "xoption": ("generator", intersectdims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison"],
    "xlabel": ["dim", "dim", "dim"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)"]
   },

   {"name": "intersect_ellipsoid",
    "executable": "benchmark_intersect",
    "config": [       
       {
          "const_configs": [],
          "fun_configs": ["Ellipsoid_intersectCoord=cached_ref", "Ellipsoid_intersectCoord=cached_reord_fma", "Ellipsoid_intersectCoord=cached_reord_fma,Ellipsoid_cacheUpdateCoord=c", "Ellipsoid_intersectCoord=cached_reord_fma,Ellipsoid_cacheUpdateCoord=fma", "Ellipsoid_intersectCoord=cached_reord_fma,Ellipsoid_cacheUpdateCoord=vec","Ellipsoid_intersectCoord=cached_reord_fma,Ellipsoid_cacheUpdateCoord=vec_u2","Ellipsoid_intersectCoord=cached_reord_fma,Ellipsoid_cacheUpdateCoord=vec_u4","Ellipsoid_intersectCoord=cached_reord_fma,Ellipsoid_cacheUpdateCoord=vec2","Ellipsoid_intersectCoord=cached_reord_fma,Ellipsoid_cacheUpdateCoord=vec2_u2","Ellipsoid_intersectCoord=cached_reord_fma,Ellipsoid_cacheUpdateCoord=vec2_u4"],
          "run_configs": ["r=100000,intersect=intersectCoord"],
          "input_configs": [("generator", intersectEbodies)]
       }
    ],
    "xoption": ("generator", intersectEdims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison"],
    "xlabel": ["dim", "dim", "dim"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)"]
   },

   {"name": "sparse_polytope_volume",
    "executable": "benchmark_A1_volume",
    "config": [ 
       {
          "const_configs": ["step_size=1000"],
          "fun_configs": [],
          "run_configs": ["r=1,polytopeType=3"],
          "input_configs": [("generator", intersectSparseBodies)]
       },
       {
          "const_configs": ["step_size=1000"],
          "fun_configs": ["PolytopeCSC_intersectCoord=cached_ref"],
          "run_configs": ["r=1,polytopeType=2"],
          "input_configs": [("generator", intersectSparseBodies)]
       },
    ],
    "xoption": ("generator", intersectSparseDims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison"],
    "xlabel": ["dim", "dim", "dim"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)"]
   },

   {"name": "sparse_polytope_intersect",
    "executable": "benchmark_intersect",
    "config": [       
       {
          "const_configs": [],
          "fun_configs": [],
          "run_configs": ["r=100000,polytopeType=3,intersect=intersectCoord"],
          "input_configs": [("generator", intersectSparseBodies)]
       },
       {
          "const_configs": [],
          "fun_configs": ["PolytopeCSC_intersectCoord=cached_ref"],#,"PolytopeCSC_intersectCoord=ref"],
          "run_configs": ["r=100000,polytopeType=2,intersect=intersectCoord"],
          "input_configs": [("generator", intersectSparseBodies)]
       },
    ],
    "xoption": ("generator", intersectSparse4Dims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison"],
    "xlabel": ["dim", "dim", "dim"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)"]
   },

   {"name": "sparse_polytopeJIT_intersect",
    "executable": "benchmark_intersect",
    "config": [       
       {
          "const_configs": [],
          "fun_configs": ["PolytopeJIT_gen=single_rax","PolytopeJIT_gen=single_data","PolytopeJIT_gen=single_data_acc","PolytopeJIT_gen=double_data","PolytopeJIT_gen=double_data,polytopeOptimize=true","PolytopeJIT_gen=quad_data,polytopeOptimize=true","PolytopeJIT_gen=quad_data","PolytopeJIT_gen=quad_data_acc,polytopeOptimize=true"],
          "run_configs": ["r=100000,polytopeType=3,intersect=intersectCoord_only"],
          "input_configs": [("generator", intersectSparse4Bodies)]
       },
    ],
    "xoption": ("generator", intersectSparse4Dims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison"],
    "xlabel": ["dim", "dim", "dim"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)"]
   },

   {"name": "dense_polytopeJIT_intersect",
    "executable": "benchmark_intersect",
    "config": [       
       {
          "const_configs": [],
          "fun_configs": ["PolytopeJIT_gen=single_rax","PolytopeJIT_gen=single_data","PolytopeJIT_gen=single_data_acc","PolytopeJIT_gen=double_data","PolytopeJIT_gen=double_data,polytopeOptimize=true","PolytopeJIT_gen=quad_data,polytopeOptimize=true","PolytopeJIT_gen=quad_data","PolytopeJIT_gen=quad_data_acc,polytopeOptimize=true"],
          "run_configs": ["r=100000,polytopeType=3,intersect=intersectCoord_only"],
          "input_configs": [("generator", cubeRotBodies)]
       },
    ],
    "xoption": ("generator", cubeRotDims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison"],
    "xlabel": ["dim", "dim", "dim"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)"]
   },

   {"name": "sparse_polytopeJIT_update",
    "executable": "benchmark_intersect",
    "config": [       
       {
          "const_configs": [],
          "fun_configs": ["PolytopeJIT_gen=single_rax","PolytopeJIT_gen=single_data","PolytopeJIT_gen=single_data_acc","PolytopeJIT_gen=double_data","PolytopeJIT_gen=double_data,polytopeOptimize=true","PolytopeJIT_gen=quad_data,polytopeOptimize=true","PolytopeJIT_gen=quad_data","PolytopeJIT_gen=quad_data_acc,polytopeOptimize=true"],
          "run_configs": ["r=200000,polytopeType=3,intersect=cacheUpdateCoord"],
          "input_configs": [("generator", intersectSparse4Bodies)]
       },
    ],
    "xoption": ("generator", intersectSparse4Dims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison"],
    "xlabel": ["dim", "dim", "dim"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)"]
   },

   {"name": "dense_polytopeJIT_update",
    "executable": "benchmark_intersect",
    "config": [       
       {
          "const_configs": [],
          "fun_configs": ["PolytopeJIT_gen=single_rax","PolytopeJIT_gen=single_data","PolytopeJIT_gen=single_data_acc","PolytopeJIT_gen=double_data","PolytopeJIT_gen=double_data,polytopeOptimize=true","PolytopeJIT_gen=quad_data,polytopeOptimize=true","PolytopeJIT_gen=quad_data","PolytopeJIT_gen=quad_data_acc,polytopeOptimize=true"],
          "run_configs": ["r=200000,polytopeType=3,intersect=cacheUpdateCoord"],
          "input_configs": [("generator", cubeRotBodies)]
       },
    ],
    "xoption": ("generator", cubeRotDims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison"],
    "xlabel": ["dim", "dim", "dim"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)"]
   },

   {"name": "jit_test",
    "executable": "benchmark_jit",
    "config": [       
       {
          "const_configs": [],
          "fun_configs": [],
          "run_configs": ["r=200000,experiment=test","r=100000,experiment=test2"],
          "input_configs": [("n", jitTest)]
       },
    ],
    "xoption": ("n", jitTestDims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison"],
    "xlabel": ["n", "n", "n"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)"]
   },

   {"name": "jit_test_3",
    "executable": "benchmark_jit",
    "config": [       
       {
          "const_configs": [],
          "fun_configs": [],
          "run_configs": [
              "r=100000,experiment=test3,m=1,w=1",
              "r=100000,experiment=test3,m=10,w=2",
              "r=100000,experiment=test3,m=100,w=3",
              "r=100000,experiment=test3,m=100,w=2",
              "r=100000,experiment=test3,m=100,w=1",
              "r=100000,experiment=test3,m=10,w=3",
              "r=100000,experiment=test3,m=1,w=3",
              ],
          "input_configs": [("n", jitTest)]
       },
    ],
    "xoption": ("n", jitTestDims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison"],
    "xlabel": ["n", "n", "n"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)"]
   },

   {"name": "jit_test_4",
    "executable": "benchmark_jit",
    "config": [       
       {
          "const_configs": [],
          "fun_configs": [],
          "run_configs": [
              "r=100000,experiment=test4,m=1,w=1",
              "r=100000,experiment=test4,m=10,w=2",
              "r=100000,experiment=test4,m=100,w=3",
              "r=100000,experiment=test4,m=100,w=2",
              "r=100000,experiment=test4,m=100,w=1",
              "r=100000,experiment=test4,m=10,w=3",
              "r=100000,experiment=test4,m=1,w=3",
              ],
          "input_configs": [("n", jitTest)]
       },
    ],
    "xoption": ("n", {k:int(v)/2 for (k,v) in jitTestDims.items()}),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison"],
    "xlabel": ["n", "n", "n"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)"]
   },


   {"name": "cube_rot_volume",
    "executable": "benchmark_A1_volume",
    "config": [ 
       {
          "const_configs": ["step_size=1000"],
          "fun_configs": [],
          "run_configs": ["r=1,polytopeType=3"],
          "input_configs": [("generator", cubeRotBodies)]
       },
       {
          "const_configs": ["step_size=1000"],
          "fun_configs": ["PolytopeCSC_intersectCoord=cached_ref", "PolytopeCSC_intersectCoord=cached_b_ref"],
          "run_configs": ["r=1,polytopeType=2"],
          "input_configs": [("generator", cubeRotBodies)]
       },
       {
          "const_configs": ["step_size=1000"],
          "fun_configs": [],
          "run_configs": ["r=1,polytopeType=4"],
          "input_configs": [("generator", cubeRotBodies)]
       },
       {
          "const_configs": ["step_size=1000"],
          "fun_configs": ["Polytope_intersectCoord=cached_ref"],
          "run_configs": ["r=1,polytopeType=0"],
          "input_configs": [("generator", cubeRotBodies)]
       },
       {
          "const_configs": ["step_size=1000"],
          "fun_configs": ["PolytopeT_intersectCoord=cached_b_ref"],
          "run_configs": ["r=1,polytopeType=1"],
          "input_configs": [("generator", cubeRotBodies)]
       },
    ],
    "xoption": ("generator", cubeRotDims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison"],
    "xlabel": ["dim", "dim", "dim"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)"]
   },
   
   {"name": "csc_intersect_9_sparse_bodies_all_cache",
    "executable": "benchmark_intersect",
    "config": [
       {
          "const_configs": [],
          "fun_configs": ["PolytopeCSC_intersectCoord=cached_b_ref",
                          #"PolytopeCSC_intersectCoord=ref",
                          "PolytopeCSC_intersectCoord=cached_b_vec",
                          "PolytopeCSC_intersectCoord=cached_b_vec_vec_nan",
                          "PolytopeCSC_intersectCoord=cached_ref",
                          "PolytopeCSC_intersectCoord=cached_b_vec_nogather",
                          "PolytopeCSC_intersectCoord=cached_b_vec_nan",
                          "PolytopeCSC_intersectCoord=cached_b_vec_nan_inv",
                          #"PolytopeCSC_intersectCoord=cached_b_vec_inl",
                          #"PolytopeCSC_intersectCoord=cached_b_vec_inl_2accs",
          ],
          "run_configs": ["r=100000,polytopeType=2,intersect=intersectCoord_only"],
          "input_configs": [("generator", intersectSparseBodies)]
       },
    ],
    "xoption": ("generator", intersectSparseDims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison"],
    "xlabel": ["dim", "dim", "dim"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)"]
   },

   
   {"name": "csc_intersect_9_dense_bodies_all_cache",
    "executable": "benchmark_intersect",
    "config": [
       {
          "const_configs": [],
          "fun_configs": ["PolytopeCSC_intersectCoord=cached_b_ref",
                          #"PolytopeCSC_intersectCoord=ref",
                          "PolytopeCSC_intersectCoord=cached_b_vec",
                          "PolytopeCSC_intersectCoord=cached_b_vec_vec_nan",
                          "PolytopeCSC_intersectCoord=cached_ref",
                          "PolytopeCSC_intersectCoord=cached_b_vec_nogather",
                          "PolytopeCSC_intersectCoord=cached_b_vec_nan",
                          "PolytopeCSC_intersectCoord=cached_b_vec_nan_inv",
                          #"PolytopeCSC_intersectCoord=cached_b_vec_inl",
                          #"PolytopeCSC_intersectCoord=cached_b_vec_inl_2accs",
          ],
          "run_configs": ["r=100000,polytopeType=2,intersect=intersectCoord_only"],
          "input_configs": [("generator", mbintersectbodies)]
       },
    ],
    "xoption": ("generator", mbintersectdims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison"],
    "xlabel": ["dim", "dim", "dim"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)"]
   },

   
   {"name": "csc_readtest_sparse",
    "executable": "benchmark_intersect",
    "config": [
       {
          "const_configs": [],
          "fun_configs": ["PolytopeCSC_intersectCoord=cached_vec_onlyread",
                          "PolytopeCSC_intersectCoord=cached_b_vec_nan_inv",
          ],
          "run_configs": ["r=100000,polytopeType=2,intersect=intersectCoord_only"],
          "input_configs": [("generator", intersectSparseBodies)]
       },
    ],
    "xoption": ("generator", intersectSparseDims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison"],
    "xlabel": ["dim", "dim", "dim"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)"]
   },

   
   {"name": "T_intersect_vecs",
    "executable": "benchmark_intersect",
    "config": [
       {
          "const_configs": [],
          "fun_configs": [#"PolytopeT_intersectCoord=ref",
                          #"PolytopeT_intersectCoord=cached_ref",
                          #"PolytopeT_intersectCoord=cached_b_ref",
                          "PolytopeT_intersectCoord=cached_b_vec",
                          "PolytopeT_intersectCoord=cached_b_vec2",
                          "PolytopeT_intersectCoord=cached_b_vec_inl",
                          #"PolytopeT_intersectCoord=cached_vectorized",                          
          ],
          "run_configs": ["r=100000,polytopeType=1,intersect=intersectCoord_only"],
          "input_configs": [("generator", mbintersectbodies)]
       },
    ],
    "xoption": ("generator", mbintersectdims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison"],
    "xlabel": ["dim", "dim", "dim"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)"]
   },

   {"name": "T_intersect_only",
    "executable": "benchmark_intersect",
    "config": [
       {
          "const_configs": [],
          "fun_configs": ["PolytopeT_intersectCoord=ref",
                          "PolytopeT_intersectCoord=cached_ref",
                          "PolytopeT_intersectCoord=cached_b_ref",
                          "PolytopeT_intersectCoord=cached_b_vec",
                          "PolytopeT_intersectCoord=cached_b_vec2",
                          #"PolytopeT_intersectCoord=cached_vectorized",                          
          ],
          "run_configs": ["r=100000,polytopeType=1,intersect=intersectCoord_only"],
          "input_configs": [("generator", mbintersectbodies)]
       },
    ],
    "xoption": ("generator", mbintersectdims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison"],
    "xlabel": ["dim", "dim", "dim"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)"]
   },

   {"name": "T_intersect_update",
    "executable": "benchmark_intersect",
    "config": [
       {
          "const_configs": [],
          "fun_configs": ["PolytopeT_intersectCoord=ref",
                          "PolytopeT_intersectCoord=cached_ref",
                          "PolytopeT_intersectCoord=cached_b_ref",
                          "PolytopeT_intersectCoord=cached_b_vec",
                          "PolytopeT_intersectCoord=cached_b_vec2",
                          #"PolytopeT_intersectCoord=cached_vectorized",                          
          ],
          "run_configs": ["r=100000,polytopeType=1,intersect=cacheUpdateCoord"],
          "input_configs": [("generator", mbintersectbodies)]
       },
    ],
    "xoption": ("generator", mbintersectdims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison"],
    "xlabel": ["dim", "dim", "dim"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)"]
   },
   
   {"name": "csc_cache_update",
    "executable": "benchmark_intersect",
    "config": [
       {
          "const_configs": [],
          "fun_configs": ["PolytopeCSC_intersectCoord=cached_b_ref",
                          "PolytopeCSC_intersectCoord=ref",
                          "PolytopeCSC_intersectCoord=cached_b_vec",
          ],
          "run_configs": ["r=100000,polytopeType=2,intersect=cacheUpdateCoord"],
          "input_configs": [("generator", intersectbodies)]
       },
    ],
    "xoption": ("generator", intersectdims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison"],
    "xlabel": ["dim", "dim", "dim"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)"]
   },

   {"name": "randomness",
    "executable": "benchmark_randomness",
    "config": [       
       {
          "const_configs": [],
          "fun_configs": ["rand_f=std_rand", "rand_f=std_rand_chunked", "rand_f=sr_rand", "rand_f=sr_rand_chunked", "rand_f=sr_rand_vec", "rand_f=mt_rand"],
          "run_configs": ["r=10000,d=16384,i=16384,rand_chunk_size=512"],
          "input_configs": [("rand_val_t", randValTypes)]
       }
    ],
    "xoption": ("rand_val_t", randValTypes_),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison"],
    "xlabel": ["rand value type", "rand value type", "rand value type"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)"]
   },

   {"name": "get_rd_0_1",
    "executable": "benchmark_randomness",
    "config": [       
       {
          "const_configs": [],
          "fun_configs": rd_0_1,
          "run_configs": ["r=10000,rand_val_t=random_double_0_1"],
          "input_configs": [("d", [2**i for i in range(10,16)])],
       }
    ],
    "xoption": ("d", {str(2**i): str(2**i) for i in range(10,16)}),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison"],
    "xlabel": ["#doubles", "#doubles", "#doubles"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)"]
   },

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
            print("     ",b["name"]);
         sys.exit(1);


         
"""
config contains the following:
const_config :: [string]
fun_config :: [string]
run_config :: [string]
input_conifg :: [(char, [(readable, char)])]

from this we create the following list of config strings
'-c "const" -f "<fun>" -n "<run,c=str(cconf)>"' 
for 
   const : const_config, 
   fun : fun_config, 
   run : run_config,
   (c, conf) : input_config
   cconf : conf
"""
def get_config(config):
   consts = [('-c', const) for const in config["const_configs"]] or [('','')]
   funs = [('-f', fun) for fun in config["fun_configs"]] or [('','')]
   inputs = [('-b', ','.join(c_prod)) for c_prod in
             itertools.product(
                *(list(map(
                      lambda c: ['{}={}'.format(c[0], conf) for conf in c[1]],
                      config["input_configs"]
                )) + [config["run_configs"]])
             )
   ] or [('','')]
   return itertools.product(consts, funs, inputs)


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
def run_benchmark(bname, bexe, config_strings):
   results = []

   for constconf,runconf,inputconf in config_strings:
      config_string = '{} "{}" {} "{}" {} "{}"'.format(constconf[0], constconf[1], runconf[0], runconf[1], inputconf[0], inputconf[1])
      config_string_printable = (
         config_string
         .replace(" ", "_")
         .replace(",", "")
      ).replace("__","")

      print('# Running Benchmark \n{} {}\n'.format(sys.path[0]+"/"+bexe, config_string));
      myenv = os.environ;
      proc = subprocess.Popen(
         [sys.path[0]+"/"+bexe,
          constconf[0],
          constconf[1],
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
            ddd = eval(line)
            assert(line[0]==123)
            results.append((config_string, ddd))
            f.write(str(ddd)+'\n')
         except:
            f.write(line.decode('utf-8'))
      f.close()

   return results


for benchmark in DO_BENCHMARKS:
   bname = benchmark["name"]
   bexe = benchmark["executable"]
   result = run_benchmark(bname,
                          bexe,
                          get_configs(benchmark)
   )
   # get x-axis labels and add them to data 
   result = list(map(lambda res: (*res, get_label(benchmark["xoption"], res[0])), result))
   pprint.pprint(result)
   plot(sys.path[0], bname, result, benchmark["xoption"][0], benchmark["title"], benchmark["xlabel"], benchmark["ylabel"])
