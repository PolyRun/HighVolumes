#!/usr/bin/env python3

import sys
import os
import subprocess
import pprint
import operator
import itertools
import re
from threading import Timer
from plot import plot

BENCH_MAX_TIME = 180	 #seconds

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

#intersectEbodies = [ "ball_r1.0_3", "ball_r1.0_10", "ball_r1.0_20", "ball_r1.0_40", "ball_r1.0_60", "ball_r1.0_100", "ball_r1.0_150", "ball_r1.0_200"]
#intersectEdims = {"ball_r1.0_10": '10', "ball_r1.0_3": '3', "ball_r1.0_20": '20', "ball_r1.0_40": '40', "ball_r1.0_60": '60', "ball_r1.0_100": '100', "ball_r1.0_150": '150', "ball_r1.0_200": '200'}

#intersectEbodies = ["2sphere_preprocessed_3","2sphere_preprocessed_10","2sphere_preprocessed_20","2sphere_preprocessed_40","2sphere_preprocessed_60","2sphere_preprocessed_100","2sphere_preprocessed_150","2sphere_preprocessed_200","2sphere_preprocessed_250","2sphere_preprocessed_300",]
#intersectEdims = {"2sphere_preprocessed_10": '10', "2sphere_preprocessed_3": '3', "2sphere_preprocessed_20": '20', "2sphere_preprocessed_40": '40', "2sphere_preprocessed_60": '60', "2sphere_preprocessed_100": '100', "2sphere_preprocessed_150": '150', "2sphere_preprocessed_200": '200', "2sphere_preprocessed_250": '250', "2sphere_preprocessed_300": '300'}

intersectEbodies = [ "ellipsoid_3", "ellipsoid_10", "ellipsoid_20", "ellipsoid_40", "ellipsoid_60", "ellipsoid_100", "ellipsoid_150", "ellipsoid_200"]
intersectEdims = {"ellipsoid_10": '10', "ellipsoid_3": '3', "ellipsoid_20": '20', "ellipsoid_40": '40', "ellipsoid_60": '60', "ellipsoid_100": '100', "ellipsoid_150": '150', "ellipsoid_200": '200'}


crossBodies = ['cross_rn_{}'.format(i) for i in range(3,14,2)]
crossDims = {'cross_rn_{}'.format(i): i for i in range(3,14,2)}
crossBodies2 = ['cross_rn_{}'.format(i) for i in [3,4]]
crossDims2 = {'cross_rn_{}'.format(i): i for i in [3,4]}

crossPDims = [4,6,8,10,12]
crossPDims = {"cross_rn_"+str(i):str(i) for i in crossPDims}
crossPBodies = [ name for (name,i) in crossPDims.items()]

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

jitTest12 = [4*i for i in range(1,20)]
jitTest12Dims = {str(i):str(4*i) for i in jitTest12}
jitTest12 = [str(i) for i in jitTest12]

randValTypes = ["random_int", "random_int_in_range", "random_double", "random_double_0_1", "random_double_normal", "random_double_in_range"]
randValTypes_ = {"random_int": '0', "random_int_in_range": '1', "random_double": '2', "random_double_0_1": '3', "random_double_normal": '4', "random_double_in_range": '5'}

# maybe add larger ones later
csc_jit_dims_2var = {'2var_{}'.format(i): str(i) for i in [4,5,10,20,40,60,100]} #,150,200]}
csc_jit_bodies_2var = [name for name in csc_jit_dims_2var.keys()]
csc_jit_dims_4var = {'4var_{}'.format(i): str(i) for i in [4,5,10,20,40,60,100]} #,150,200]}
csc_jit_bodies_4var = [name for name in csc_jit_dims_4var.keys()]
csc_jit_dims_2var_pre = {'2var_TSP_{}'.format(i): str(i) for i in [4,5,10,20,30,40,50,60,100]} #,150,200]}
csc_jit_bodies_2var_pre = [name for name in csc_jit_dims_2var_pre.keys()]
csc_jit_dims_cube_rot = {'cube_rot_r1.0_{}'.format(i): str(i) for i in [3,10,20,40,60,80,100]}
csc_jit_bodies_cube_rot = [name for name in csc_jit_dims_cube_rot.keys()]
csc_jit_dims_cross_rot = {'cross_rot_rn_{}'.format(i): str(i) for i in range(3,14,2)}
csc_jit_bodies_cross_rot = [name for name in csc_jit_dims_cross_rot.keys()]

csc_jit_dims_2var_vinci = {'vinci_2var_TSP_{}.ine'.format(i): str(i) for i in [10,20,30,40,50,60,100,150,200]}
csc_jit_bodies_2var_vinci = [name for name in csc_jit_dims_2var_vinci.keys()]

csc_jit_dims_dens = {'dens_{}'.format(i): str(min(int(1.55 ** i),100)) for i in range(1,12)}
csc_jit_bodies_dens = [name for name in csc_jit_dims_dens.keys()]
csc_jit_dims_dens200 = {'dens200_{}'.format(i): str(min(int(1.55 ** i),100)) for i in range(1,12)}
csc_jit_bodies_dens200 = [name for name in csc_jit_dims_dens200.keys()]

csc_jit_bodies = [
   ("2var", csc_jit_bodies_2var, csc_jit_dims_2var),
   ("4var", csc_jit_bodies_4var, csc_jit_dims_4var),
   ("2var_vinci", csc_jit_bodies_2var_vinci, csc_jit_dims_2var_vinci),
   ("cube_rot", csc_jit_bodies_cube_rot, csc_jit_dims_cube_rot),
   ("cross", csc_jit_bodies_cross_rot, csc_jit_dims_cross_rot),
   ("dens", csc_jit_bodies_dens, csc_jit_dims_dens),
]



intersects_funs = [
   ("cacheUpdateCoord",
    [
       (
          ['PolytopeCSC_intersectCoord={}'.format(fun) for fun in
           [
              "cached_b_ref",
              #"cached_b_vec",
              "cached_b_vec_vec_nan"
           ]
          ],
          2
       ),
       (
          ['PolytopeJIT_gen={}'.format(fun) for fun in
           [
              #"single_rax",
              "single_data",
              "double_data",
              "quad_data"
           ]       
          ],
          3
       )
    ]
   ),
   ("intersectCoord_only",
    [
       (
          ['PolytopeCSC_intersectCoord={}'.format(fun) for fun in
           [
              #"cached_b_ref",
              "cached_b_vec",
              "cached_b_vec_nogather",
              "cached_b_vec_nan_inv",
              #"cached_vec_onlyread"
           ]
          ],
          2
       ),
       (
          ['PolytopeJIT_gen={}'.format(fun) for fun in
           [
              #"single_rax",
              "single_data",
              "double_data",
              "quad_data",
              #"quad_data_acc"
           ]
          ],
          3
       )
    ]
   )
]
   # - Detail analysis CSC, JIT. Diffent levels of sparsity 2var, 4var. 2var preprocessed? cross, cubeRot. - take sizes where you see performance decreasing.
   #   - separate intersect_only, cacheUpdateCoord. JIT will be better because different load instructions/patterns.


csc_jit_bm = [
   {"name": 'csc_jit_only_cacheb_{}_{}'.format(name,intersect),
    "executable": "benchmark_intersect",
    "config": [
       {
          "const_configs": [],
          "fun_configs": funs,
          "run_configs": ['r=100000,polytopeType={},intersect={},polytopeOptimize=true'.format(bodytype,intersect)],
          "input_configs": [("generator", bodies)]
       }
       for funs, bodytype in fconf
    ],
    "xoption": ("generator", dims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
   }
   for name, bodies, dims in csc_jit_bodies for intersect,fconf in intersects_funs
]



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
		 perf_roofs:    Performance roofs that should be included int plots w/o machine roof
		 mem_roofs:    Memory roofs that should be included int plots w/o machine roof
		 
'''
BENCHMARKS = [

   
   
   # Listing what to plot:
   # - On different bodies compare bodies. For intersect with update. And for volume (incl polyvest) - random versions?. Only best functions.
   #   -> also roofline.
   #   -> Runtime: PolytopeT (base-cached-b, vectorized, [inv], random, [squaredNorm-cached]) vs CSC (optimized) vs JIT (optimized) vs Polyvest
   # - Detail analysis CSC, JIT. Diffent levels of sparsity 2var, 4var. 2var preprocessed? cross, cubeRot. - take sizes where you see performance decreasing.
   #   - separate intersect_only, cacheUpdateCoord. JIT will be better because different load instructions/patterns.
   # - CSC/JIT benchmark-tests? CSC extra function. JIT tests??? JIT test1/2 -> loadsd to get roof.

   # TODO:
   # choose best of PolytopeT, CSC, JIT, all cached with b.
   # Check crossBodies dimensions
   # do same for sparse bodies (2var)
   # do same for cubeRot
   # all of these also with cacheUpdateCoord

   # TODO: Jonathan: fix/use/ or delete
   {"name": "polytopeT_intersectCoord",
    "executable": "benchmark_intersect",
    "config": [       
       {
          "const_configs": [],
          #"fun_configs": ["PolytopeT_intersectCoord=cached_nc1", "PolytopeT_intersectCoord=ref", "PolytopeT_intersectCoord=cached_b_vec"],
          "fun_configs": [
                          "PolytopeT_intersectCoord=ref",
                          "PolytopeT_intersectCoord=cached_b_ref",
                          "PolytopeT_intersectCoord=cached_b_inv_vec",
                          ],
          "run_configs": ["intersect=intersectCoord,polytopeType=1,r=100000"],
          "input_configs": [("generator", mbintersectbodies)]
       },
    ],
    "xoption": ("generator", mbintersectdims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [20.0/27, 5.0],
    "mem_roofs": []
   },

   {"name": "polytopeT_intersectCoord_nobase",
    "executable": "benchmark_intersect",
    "config": [       
       {
          "const_configs": [],
          #"fun_configs": ["PolytopeT_intersectCoord=cached_nc1", "PolytopeT_intersectCoord=ref", "PolytopeT_intersectCoord=cached_b_vec"],
          "fun_configs": [
                          #"PolytopeT_intersectCoord=ref",
                          "PolytopeT_intersectCoord=cached_ref",
                          "PolytopeT_intersectCoord=cached_b_ref",
                          "PolytopeT_intersectCoord=cached_b_vec",
                          #"PolytopeT_intersectCoord=cached_b_vec2",
                          #"PolytopeT_intersectCoord=cached_b_vec_inl"
                          "PolytopeT_intersectCoord=cached_b_inv_ref",
                          "PolytopeT_intersectCoord=cached_b_inv_vec",
                          ],
          "run_configs": ["intersect=intersectCoord_only,polytopeType=1,r=100000"],
          "input_configs": [("generator", mbintersectbodies)]
       },
    ],
    "xoption": ("generator", mbintersectdims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [20.0/27],
    "mem_roofs": []
   },

   {"name": "polytopeT_update",
    "executable": "benchmark_intersect",
    "config": [       
       {
          "const_configs": [],
          #"fun_configs": ["PolytopeT_intersectCoord=cached_nc1", "PolytopeT_intersectCoord=ref", "PolytopeT_intersectCoord=cached_b_vec"],
          "fun_configs": [
                          "PolytopeT_intersectCoord=ref",
                          "PolytopeT_intersectCoord=cached_ref",
                          "PolytopeT_intersectCoord=cached_b_ref",
                          "PolytopeT_intersectCoord=cached_b_vec",
                          #"PolytopeT_intersectCoord=cached_b_vec2",
                          #"PolytopeT_intersectCoord=cached_b_vec_inl"
                          "PolytopeT_intersectCoord=cached_b_inv_ref",
                          "PolytopeT_intersectCoord=cached_b_inv_vec",
                          ],
          "run_configs": ["intersect=cacheUpdateCoord,polytopeType=1,r=100000"],
          "input_configs": [("generator", mbintersectbodies)]
       },
    ],
    "xoption": ("generator", mbintersectdims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [20.0/27],
    "mem_roofs": []
   },

   {"name": "polytopeT_update_nobase",
    "executable": "benchmark_intersect",
    "config": [       
       {
          "const_configs": [],
          #"fun_configs": ["PolytopeT_intersectCoord=cached_nc1", "PolytopeT_intersectCoord=ref", "PolytopeT_intersectCoord=cached_b_vec"],
          "fun_configs": [
                          #"PolytopeT_intersectCoord=ref",
                          "PolytopeT_intersectCoord=cached_ref",
                          "PolytopeT_intersectCoord=cached_b_ref",
                          "PolytopeT_intersectCoord=cached_b_vec",
                          #"PolytopeT_intersectCoord=cached_b_vec2",
                          #"PolytopeT_intersectCoord=cached_b_vec_inl"
                          "PolytopeT_intersectCoord=cached_b_inv_ref",
                          "PolytopeT_intersectCoord=cached_b_inv_vec",
                          ],
          "run_configs": ["intersect=cacheUpdateCoord,polytopeType=1,r=100000"],
          "input_configs": [("generator", mbintersectbodies)]
       },
    ],
    "xoption": ("generator", mbintersectdims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [20.0/27],
    "mem_roofs": []
   },


   # TODO later
   
   {"name": "ellipsoid_intersectCoord_only",
    "executable": "benchmark_intersect",
    "config": [       
       {
          "const_configs": [],
          "fun_configs": ["Ellipsoid_intersectCoord=cached_ref","Ellipsoid_intersectCoord=cached_reord","Ellipsoid_intersectCoord=cached_reord2","Ellipsoid_intersectCoord=cached_reord3","Ellipsoid_intersectCoord=cached_reord_fma"],
          "run_configs": ["r=100000,intersect=intersectCoord_only"],
          "input_configs": [("generator", intersectEbodies)]
       }
    ],
    "xoption": ("generator", intersectEdims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
   },
   
   {"name": "ellipsoid_cacheUpdateCoord",
    "executable": "benchmark_intersect",
    "config": [       
       {
          "const_configs": [],
          "fun_configs": ["Ellipsoid_cacheUpdateCoord=ref","Ellipsoid_cacheUpdateCoord=fma","Ellipsoid_cacheUpdateCoord=vec","Ellipsoid_cacheUpdateCoord=vec_u2"],
          "run_configs": ["r=100000,intersect=cacheUpdateCoord"],
          "input_configs": [("generator", intersectEbodies)]
       }
    ],
    "xoption": ("generator", intersectEdims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
   },
   
   {"name": "ellipsoid_intersect",
    "executable": "benchmark_intersect",
    "config": [       
       {
          "const_configs": [],
          "fun_configs": ["Ellipsoid_cacheUpdateCoord=ref,Ellipsoid_intersectCoord=cached_ref","Ellipsoid_cacheUpdateCoord=fma,Ellipsoid_intersectCoord=cached_reord_fma","Ellipsoid_cacheUpdateCoord=vec_u2,Ellipsoid_intersectCoord=cached_reord_fma"],
          "run_configs": ["r=100000,intersect=intersectCoord"],
          "input_configs": [("generator", intersectEbodies)]
       }
    ],
    "xoption": ("generator", intersectEdims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
   },
   
   {"name": "ellipsoid_volume",
    "executable": "benchmark_A1_volume",
    "config": [       
       {
          "const_configs": [],
          "fun_configs": ["Ellipsoid_cacheUpdateCoord=ref,Ellipsoid_intersectCoord=cached_ref","Ellipsoid_cacheUpdateCoord=fma,Ellipsoid_intersectCoord=cached_reord_fma","Ellipsoid_cacheUpdateCoord=vec_u2,Ellipsoid_intersectCoord=cached_reord_fma"],
          "run_configs": ["r=1,intersect=intersectCoord"],
          "input_configs": [("generator", intersectEbodies)]
       }
    ],
    "xoption": ("generator", intersectEdims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
   },
   
   # determine basic read roof line
   {"name": "jit_test12",
    "executable": "benchmark_jit",
    "config": [       
       {
          "const_configs": [],
          "fun_configs": [],
          "run_configs": ["r=200000,experiment=test","r=100000,experiment=test2","r=200000,experiment=testa","r=100000,experiment=test2a"],
          "input_configs": [("n", jitTest12)]
       },
    ],
    "xoption": ("n", jitTest12Dims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["n", "n", "n", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
   },
   
   # determine basic latency for kernel functions
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
              "r=100000,experiment=test3,m=200,w=1",
              "r=100000,experiment=test3,m=10,w=3",
              "r=100000,experiment=test3,m=1,w=3",
              ],
          "input_configs": [("n", jitTest)]
       },
    ],
    "xoption": ("n", jitTestDims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["n", "n", "n", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
   },

   # TODO use later?
   {"name": "randomness",
    "executable": "benchmark_randomness",
    "config": [       
       {
          "const_configs": [],
          "fun_configs": ["rand_f=std_rand", "rand_f=std_rand_chunked", "rand_f=sr_rand", "rand_f=sr_rand_chunked", "rand_f=sr_rand_vec", "rand_f=mt_rand"],
          "run_configs": ["r=10000,rand_val_t=random_int"],
          "input_configs": [("i", [2**i for i in range(14,18)])]
       }
    ],
    "xoption": ("i", {str(2**i): str(2**i) for i in range(14,18)}),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["#random ints", "#random ints", "#random ints", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
   },

]

# - On different bodies compare bodies. For intersect with update. And for volume (incl polyvest) - random versions?. Only best functions.
#   -> also roofline.
#   -> Runtime: PolytopeT (base-cached-b, vectorized, [inv], random, [squaredNorm-cached]) vs CSC (optimized) vs JIT (optimized) vs Polyvest

def pushConfig(func,bodyname,bodies,dims,r,step_size):
   executable = {
           "volume":"benchmark_A1_volume",
           "intersect":"benchmark_intersect",
           "update":"benchmark_intersect",
           }[func]
   runConf = {
           "volume":"",
           "intersect":",intersect=intersectCoord_only",
           "update":",intersect=cacheUpdateCoord",
           }[func]

   config = [ 
       #{
       #   "const_configs": ["step_size="+str(step_size)],
       #   "fun_configs": [
       #       "PolytopeT_intersectCoord=cached_b_ref"
       #       ],
       #   "run_configs": ["r="+str(r)+",polytopeType=1"+runConf],
       #   "input_configs": [("generator", bodies)]
       #},
       #{
       #   "const_configs": ["step_size="+str(step_size)],
       #   "fun_configs": [
       #       "PolytopeT_intersectCoord=cached_b_vec"
       #       ],
       #   "run_configs": ["r="+str(r)+",polytopeType=1"+runConf],
       #   "input_configs": [("generator", bodies)]
       #},
       #{
       #   "const_configs": ["step_size="+str(step_size)],
       #   "fun_configs": [
       #       "PolytopeT_intersectCoord=cached_b_inv_ref"
       #       ],
       #   "run_configs": ["r="+str(r)+",polytopeType=1"+runConf],
       #   "input_configs": [("generator", bodies)]
       #},
       {
          "const_configs": ["step_size="+str(step_size)],
          "fun_configs": [
              "PolytopeT_intersectCoord=cached_b_inv_vec"
              ],
          "run_configs": ["r="+str(r)+",polytopeType=1"+runConf],
          "input_configs": [("generator", bodies)]
       },
       {
          "const_configs": ["step_size="+str(step_size)],
          "fun_configs": [
              "PolytopeT_intersectCoord=cached_b_inv_vec,walk_f=walkCoord_{}".format(i) for i in [1,4,8] # may lead to issues for things other than volume
              ],
          "run_configs": ["r="+str(r)+",polytopeType=1"+runConf],
          "input_configs": [("generator", bodies)]
       },
       #{
       #   "const_configs": ["step_size="+str(step_size)],
       #   "fun_configs": [
       #       "PolytopeT_intersectCoord=cached_b_vec,rand_f=std_rand"
       #       ],
       #   "run_configs": ["r="+str(r)+",polytopeType=1"+runConf],
       #   "input_configs": [("generator", bodies)]
       #},
       {
          "const_configs": ["step_size="+str(step_size)],
          "fun_configs": [],
          "fun_configs": ["PolytopeCSC_intersectCoord=cached_b_vec_nan_inv"],
          "run_configs": ["r=100,polytopeType=2"],
          "run_configs": ["r="+str(r)+",polytopeType=2"+runConf],
          "input_configs": [("generator", bodies)]
       },
       {
          "const_configs": ["step_size="+str(step_size)],
          "fun_configs": [],
          "fun_configs": ["PolytopeJIT_intersectCoord=quad_data"], # acc?
          "run_configs": ["r="+str(r)+",polytopeType=3"+runConf],
          "input_configs": [("generator", bodies)]
       }
       ]
    
   if(func == "volume"): # add Polyvest
      config.append({
          "const_configs": ["step_size="+str(step_size)],
          "fun_configs": [],
          "run_configs": ["r="+str(r)+",polytopeType=4"+runConf],
          "input_configs": [("generator", bodies)]
      }) 

   BENCHMARKS.append(
   {"name": "polytopes_" + func + "_" + bodyname,
    "executable": executable,
    "config": config,
    "xoption": ("generator", dims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": [],    
   })

## Actually add the comparison configs:
rFactor = 1
sFactor = 1
pushConfig("intersect","cubeRot",cubeRotBodies,cubeRotDims,10000*rFactor,1600)
pushConfig("intersect","2var",intersectSparseBodies,intersectSparseDims,3000*rFactor,1600)
pushConfig("intersect","4var",intersectSparse4Bodies,intersectSparse4Dims,2000*rFactor,1600)
pushConfig("intersect","cross",crossPBodies,crossPDims,10000*rFactor,1600)
#pushConfig("intersect","",,,100)
pushConfig("update","cubeRot",cubeRotBodies,cubeRotDims,10000*rFactor,1600)
pushConfig("update","2var",intersectSparseBodies,intersectSparseDims,3000*rFactor,1600)
pushConfig("update","4var",intersectSparse4Bodies,intersectSparse4Dims,2000*rFactor,1600)
pushConfig("update","cross",crossPBodies,crossPDims,10000*rFactor,1600)

pushConfig("volume","cubeRot",cubeRotBodies,cubeRotDims,1,400*sFactor)
pushConfig("volume","2var",intersectSparseBodies,intersectSparseDims,1,400*sFactor)
pushConfig("volume","4var",intersectSparse4Bodies,intersectSparse4Dims,1,400*sFactor)
pushConfig("volume","cross",crossPBodies,crossPDims,1,400*sFactor)


# polytopeT development plot (all body types?) cross_13 performance breaks down?! let's hope
# - PolytopeT ref implementation
# - cached impl
# - cached with b
# - cached withb inverse
# - cached with b inverse and vectorized

# sparsity plots:
# 1) super sparse body (2var) -> runtime comparison with polytopeT best version, show speedup

# 2) semi-sparse (2var_TSP?) -> csc, jit, polytopeT, show breakdown of JIT
# - large semi-sparse: first jit best, then at certain dim performance degration and csc scales better and polytopeT is slower



polytopeT_bm = [
   {"name": 'polytopeT_{}_bm'.format(name),
    "executable": "benchmark_intersect",
    "config": [       
       {
          "const_configs": [],
          "fun_configs": [
                          "PolytopeT_intersectCoord=cached_b_ref",
                          "PolytopeT_intersectCoord=cached_b_vec",
                          "PolytopeT_intersectCoord=cached_b_inv_ref",
                          "PolytopeT_intersectCoord=cached_b_inv_vec",
                          ],
          "run_configs": ["intersect=intersectCoord_only,polytopeType=1,r=100000"],
          "input_configs": [("generator", bodies)]
       },       
       {
          "const_configs": [],
          "fun_configs": ["PolytopeT_intersectCoord=cached_b_inv_vec,walk_f=walkCoord_4"],
          "run_configs": ["intersect=intersectCoord_only,polytopeType=1,r=100000,intersectSet=4"],
          "input_configs": [("generator", bodies)]
       },       
       {
          "const_configs": [],
          "fun_configs": ["PolytopeT_intersectCoord=cached_b_inv_vec,walk_f=walkCoord_8"],
          "run_configs": ["intersect=intersectCoord_only,polytopeType=1,r=100000,intersectSet=8"],
          "input_configs": [("generator", bodies)]
       },
    ],
    "xoption": ("generator", dims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [20.0/27, 5.0],
    "mem_roofs": []
   }
   for bodies,dims,name in [
         (intersectbodies, intersectdims, "cubeRot"),
         (crossBodies, crossDims, "cross"),
   ]
]

# density benchmark
intersects_funs = [
   (
      ["PolytopeT_intersectCoord=cached_b_inv_vec,walk_f=walkCoord_1"],
      1,
      ["intersectSet=1"]
   ),
   (
      ['PolytopeCSC_intersectCoord={}'.format(fun) for fun in
       [
          "cached_b_ref",
          "cached_b_vec_nogather",
       ]
      ],
      2,
      [""]
   ),
   (
      ['PolytopeCSC_intersectCoord=cached_b_vec_nan_inv'],
      2,
      ["intersectSet={}".format(i) for i in [1]]
   ),
   (
      ['PolytopeJIT_gen={}'.format(fun) for fun in
       [
          "single_data",
          "quad_data",
       ]
      ],
      3,
      ["intersectSet=1"]
   )
]
# - Detail analysis CSC, JIT. Diffent levels of sparsity 2var, 4var. 2var preprocessed? cross, cubeRot. - take sizes where you see performance decreasing.
#   - separate intersect_only, cacheUpdateCoord. JIT will be better because different load instructions/patterns.


density_runtime_bm = [
   {"name": 'density_runtime_bm',
    "executable": "benchmark_intersect",
    "config": [
       {
          "const_configs": [],
          "fun_configs": funs,
          "run_configs": ['r=2000,polytopeType={},intersect=intersectCoord_only,polytopeOptimize=true,{}'.format(bodytype,bconf)],
          "input_configs": [("generator", csc_jit_bodies_dens200)]
       }
       for funs, bodytype, bconfs in intersects_funs for bconf in bconfs
    ],
    "xoption": ("generator", csc_jit_dims_dens200),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
   }
]

# density8 benchmark
intersects_funs8 = [
   (
      ["PolytopeT_intersectCoord=cached_b_inv_vec,walk_f=walkCoord_1"],
      1,
      ["intersectSet=8"]
   ),
   (
      ['PolytopeCSC_intersectCoord=cached_b_vec_nan_inv'],
      2,
      ["intersectSet=8"]
   ),
   (
      ['PolytopeJIT_gen={}'.format(fun) for fun in
       [
          "quad_data",
       ]
      ],
      3,
      ["intersectSet=8"]
   )
]
# - Detail analysis CSC, JIT. Diffent levels of sparsity 2var, 4var. 2var preprocessed? cross, cubeRot. - take sizes where you see performance decreasing.
#   - separate intersect_only, cacheUpdateCoord. JIT will be better because different load instructions/patterns.


density8_runtime_bm = [
   {"name": 'density8_runtime_bm',
    "executable": "benchmark_intersect",
    "config": [
       {
          "const_configs": [],
          "fun_configs": funs,
          "run_configs": ['r=2000,polytopeType={},intersect=intersectCoord_only,polytopeOptimize=true,{}'.format(bodytype,bconf)],
          "input_configs": [("generator", csc_jit_bodies_dens200)]
       }
       for funs, bodytype, bconfs in intersects_funs8 for bconf in bconfs
    ],
    "xoption": ("generator", csc_jit_dims_dens200),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
   }
]




sparse_intersect_funs = [
   (#default
      ['PolytopeCSC_intersectCoord=cached_b_vec_nan_inv'],
      2,
      [""]
   ),
   (#cached squareNorm
      ['PolytopeCSC_intersectCoord=cached_b_vec_nan_inv,walk_f=walkCoord_1'],
      2,
      [""]
   ),
   (# cached squareNorm and fast random
      ['PolytopeCSC_intersectCoord=cached_b_vec_nan_inv,rand_f=sr_rand,rd_0_1=fast,walk_f=walkCoord_1'],
      2,
      [""]
   ),
   (#default
      ['PolytopeJIT_gen=double_data'],
      3,
      [""]
   ),
   (#cached squareNorm
      ['PolytopeJIT_gen=double_data,walk_f=walkCoord_1'],
      3,
      [""]
   ),
   (#cached squareNorm and fast random
      ['PolytopeJIT_gen=double_data,rand_f=sr_rand,rd_0_1=fast,walk_f=walkCoord_1'],
      3,
      [""]
   ),
]

runtime_2var_bm = [
   {"name": '2var_runtime_bm',
    "executable": "benchmark_A1_volume",
    "config": [
       {
          "const_configs": ["step_size=100"],
          "fun_configs": funs,
          "run_configs": ['warmup=1,r=1,polytopeType={},polytopeOptimize=true,{}'.format(bodytype,bconf)],
          "input_configs": [("generator", csc_jit_bodies_2var)]
       }
       for funs, bodytype, bconfs in sparse_intersect_funs for bconf in bconfs
    ],
    "xoption": ("generator", csc_jit_dims_2var),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
   }   
]



log_funs = [
   (
      ['PolytopeCSC_intersectCoord=cached_b_vec_nan_inv,rand_f=sr_rand,rd_0_1=fast,walk_f=walkCoord_1,log_f={}'.format(fun)
       for fun in ['nocache', 'ref', 'binary']       
      ],
      2,
      [""]
   ),
   (#cached squareNorm and fast random
      ['PolytopeJIT_gen=double_data,rand_f=sr_rand,rd_0_1=fast,walk_f=walkCoord_1,log_f={}'.format(fun)
       for fun in ['nocache', 'ref', 'binary']
      ],
      3,
      [""]
   ),
]

log_runtime_bm = [
   {"name": 'log_runtime_bm',
    "executable": "benchmark_A1_volume",
    "config": [
       {
          "const_configs": ["step_size=100"],
          "fun_configs": funs,
          "run_configs": ['warmup=1,r=1,polytopeType={},polytopeOptimize=true,{}'.format(bodytype,bconf)],
          "input_configs": [("generator", intersectSparseBodies)]
       }
       for funs, bodytype, bconfs in log_funs for bconf in bconfs
    ],
    "xoption": ("generator", intersectSparseDims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
   }   
]




# density benchmark
allbest_funs = [
   (
      ["PolytopeT_intersectCoord=cached_b_inv_vec,walk_f=walkCoord_{}".format(i) for i in [1,4,8]],
      1,
      [""]
   ),
   (
      ['PolytopeCSC_intersectCoord=cached_b_vec_nan_inv,walk_f=walkCoord_{}'.format(i) for i in [1,4,8]],
      2,
      [""]
   ),
   (
      ['PolytopeJIT_gen=quad_data,walk_f=walkCoord_{}'.format(i) for i in [1,4,8]],
      3,
      [""]
   ),
#   (
#      ['PolytopeJIT_gen={},walk_f=walkCoord_1'.format(fun) for fun in
#       [
#          "single_data",
#          "double_data",
#          "quad_data",
#       ]
#      ],
#      3,
#      [""]
#   )
]
# - Detail analysis CSC, JIT. Diffent levels of sparsity 2var, 4var. 2var preprocessed? cross, cubeRot. - take sizes where you see performance decreasing.
#   - separate intersect_only, cacheUpdateCoord. JIT will be better because different load instructions/patterns.




allbest_bodies = [
   ("cube_rot", csc_jit_bodies_cube_rot, csc_jit_dims_cube_rot),
   ("cross", csc_jit_bodies_cross_rot, csc_jit_dims_cross_rot)
]

allbest_roofline_bm = [
   {"name": 'allbest_roofline_{}_bm'.format(name),
    "executable": "benchmark_A1_volume",
    "config": [
       {
          "const_configs": ["step_size=100"],
          "fun_configs": funs,
          "run_configs": ['warmup=1,r=10,polytopeType={},polytopeOptimize=true,{}'.format(bodytype,bconf)],
          "input_configs": [("generator", body)]
       }
       for funs, bodytype, bconfs in allbest_funs for bconf in bconfs
    ],
    "xoption": ("generator", dim),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [5.0],
    "mem_roofs": []
   }
   for name, body, dim in allbest_bodies
]


BENCHMARKS += allbest_roofline_bm
BENCHMARKS += density_runtime_bm
BENCHMARKS += density8_runtime_bm
BENCHMARKS += runtime_2var_bm
BENCHMARKS += polytopeT_bm
BENCHMARKS += csc_jit_bm
BENCHMARKS += log_runtime_bm


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
def run_benchmark(bname, bexe, config_strings, xoption):
   results = []
   failed = set()
   global is_timeout

   for constconf,runconf,inputconf in config_strings:
      r = re.compile('{}=[^,]*'.format(xoption[0]))
      is_timeout = False

      
      config_string = '{} "{}" {} "{}" {} "{}"'.format(constconf[0], constconf[1], runconf[0], runconf[1], inputconf[0], inputconf[1])
      
      if ''.join(r.split(config_string)) in failed:
         continue
      
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

      def timeOut():
         global is_timeout
         print(is_timeout)
         is_timeout = True
         print("timeout - kill process!")
         proc.kill()
         outloc = proc.stdout.read().split(b'\n')
         pprint.pprint(outloc)
         
      timer = Timer(BENCH_MAX_TIME, timeOut)

      try:
         timer.start()
         proc.wait()
      except KeyboardInterrupt:
         print("user interrupt!")
         sys.exit(1)
      finally:
         timer.cancel()

        
      if is_timeout:
         failed.add(''.join(r.split(config_string)))
         print("failed string")
         pprint.pprint(failed)
         #results.append((config_string, None))
      else:
         print("we are here")
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
                          get_configs(benchmark),
                          benchmark["xoption"]
   )
   # get x-axis labels and add them to data 
   result = list(map(lambda res: (*res, get_label(benchmark["xoption"], res[0])), result))
   pprint.pprint(result)
   plot(sys.path[0], bname, result, benchmark["xoption"][0], benchmark["title"], benchmark["xlabel"], benchmark["ylabel"], benchmark["perf_roofs"], benchmark["mem_roofs"])
