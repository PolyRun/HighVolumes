

OLDBENCHMARKS = [
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
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["n", "n", "n", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
   },

    {"name": "intersectB_polytope_cross",
    "executable": "benchmark_intersect",
    "config": [       
       {
          "const_configs": [],
          "fun_configs": ["PolytopeT_intersectCoord=cached_b_vec"],
          "run_configs": ["intersect=intersectCoord_only,polytopeType=1,r=10000000","intersect=intersectCoord_only,polytopeType=1"],
          "input_configs": [("generator", crossBodies2)]
       },
    ],
    "xoption": ("generator", crossDims2),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
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
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
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
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
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
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
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
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
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
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
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
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
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
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["n", "n", "n", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
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
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
   },

   {"name": "test_eps_2",
    "executable": "benchmark_intersect",
    "config": [
       {
          "const_configs": [],
          "fun_configs": [#"PolytopeT_intersectCoord=ref",
                          #"PolytopeT_intersectCoord=cached_ref",
                          #"PolytopeT_intersectCoord=cached_b_ref",
                          #"PolytopeT_intersectCoord=cached_b_vec",
                          #"PolytopeT_intersectCoord=cached_b_vec2",
                          "PolytopeT_intersectCoord=cached_b_vec_inl",
                          #"PolytopeT_intersectCoord=cached_vectorized",                          
          ],
          "run_configs": ["r=1,polytopeType=1,intersect=intersectCoord_only"],
          "input_configs": [("generator", mbintersectbodies)]
       },
    ],
    "xoption": ("generator", mbintersectdims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
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
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
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
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
   },

   ###################################################################################################
   # CSC start
   ###################################################################################################
   {"name": "csc_cache_update_cube_rot",
    "executable": "benchmark_intersect",
    "config": [
       {
          "const_configs": [],
          "fun_configs": ["PolytopeCSC_intersectCoord=cached_b_ref", # withb, withb
                          "PolytopeCSC_intersectCoord=cached_b_vec", # fma, fma
                          "PolytopeCSC_intersectCoord=cached_b_vec_vec_nan", # vec, vec
          ],
          "run_configs": ["r=100000,polytopeType=2,intersect=cacheUpdateCoord"],
          "input_configs": [("generator", intersectbodies)]
       },
    ],
    "xoption": ("generator", intersectdims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
   },
   
   {"name": "csc_cache_update_2var",
    "executable": "benchmark_intersect",
    "config": [
       {
          "const_configs": [],
          "fun_configs": ["PolytopeCSC_intersectCoord=cached_b_ref", # withb, withb
                          "PolytopeCSC_intersectCoord=cached_b_vec", # fma, fma
                          "PolytopeCSC_intersectCoord=cached_b_vec_vec_nan", # vec, vec
          ],
          "run_configs": ["r=100000,polytopeType=2,intersect=cacheUpdateCoord"],
          "input_configs": [("generator", intersectSparseBodies)]
       },
    ],
    "xoption": ("generator", intersectSparseDims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
   },

   {"name": "csc_cache_update_4var",
    "executable": "benchmark_intersect",
    "config": [
       {
          "const_configs": [],
          "fun_configs": ["PolytopeCSC_intersectCoord=cached_b_ref", # withb, withb
                          "PolytopeCSC_intersectCoord=cached_b_vec", # fma, fma
                          "PolytopeCSC_intersectCoord=cached_b_vec_vec_nan", # vec, vec
          ],
          "run_configs": ["r=100000,polytopeType=2,intersect=cacheUpdateCoord"],
          "input_configs": [("generator", intersectSparse4Bodies)]
       },
    ],
    "xoption": ("generator", intersectSparse4Dims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
   },

   {"name": "csc_cache_update_cross",
    "executable": "benchmark_intersect",
    "config": [
       {
          "const_configs": [],
          "fun_configs": ["PolytopeCSC_intersectCoord=cached_b_ref", # withb, withb
                          "PolytopeCSC_intersectCoord=cached_b_vec", # fma, fma
                          "PolytopeCSC_intersectCoord=cached_b_vec_vec_nan", # vec, vec
          ],
          "run_configs": ["r=100000,polytopeType=2,intersect=cacheUpdateCoord"],
          "input_configs": [("generator", crossBodies)]
       },
    ],
    "xoption": ("generator", crossDims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
   },
      
   {"name": "csc_intersect_only_2var",
    "executable": "benchmark_intersect",
    "config": [
       {
          "const_configs": [],
          "fun_configs": ["PolytopeCSC_intersectCoord=cached_b_ref",
                          "PolytopeCSC_intersectCoord=ref",
                          "PolytopeCSC_intersectCoord=cached_b_vec",
                          "PolytopeCSC_intersectCoord=cached_b_vec_vec_nan",
                          "PolytopeCSC_intersectCoord=cached_ref",
                          "PolytopeCSC_intersectCoord=cached_b_vec_nogather",
                          "PolytopeCSC_intersectCoord=cached_b_vec_nan",
                          "PolytopeCSC_intersectCoord=cached_b_vec_nan_inv",
                          "PolytopeCSC_intersectCoord=cached_b_vec_inl",
                          "PolytopeCSC_intersectCoord=cached_b_vec_inl_2accs",
          ],
          "run_configs": ["r=100000,polytopeType=2,intersect=intersectCoord_only"],
          "input_configs": [("generator", intersectSparseBodies)]
       },
    ],
    "xoption": ("generator", intersectSparseDims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
   },

   {"name": "csc_intersect_only_cross",
    "executable": "benchmark_intersect",
    "config": [
       {
          "const_configs": [],
          "fun_configs": ["PolytopeCSC_intersectCoord=cached_b_ref",
                          "PolytopeCSC_intersectCoord=ref",
                          "PolytopeCSC_intersectCoord=cached_b_vec",
                          "PolytopeCSC_intersectCoord=cached_b_vec_vec_nan",
                          "PolytopeCSC_intersectCoord=cached_ref",
                          "PolytopeCSC_intersectCoord=cached_b_vec_nogather",
                          "PolytopeCSC_intersectCoord=cached_b_vec_nan",
                          "PolytopeCSC_intersectCoord=cached_b_vec_nan_inv",
                          "PolytopeCSC_intersectCoord=cached_b_vec_inl",
                          "PolytopeCSC_intersectCoord=cached_b_vec_inl_2accs",
          ],
          "run_configs": ["r=100000,polytopeType=2,intersect=intersectCoord_only"],
          "input_configs": [("generator", crossBodies)]
       },
    ],
    "xoption": ("generator", crossDims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
   },

   {"name": "csc_intersect_only_2var_cacheb_new",
    "executable": "benchmark_intersect",
    "config": [
       {
          "const_configs": [],
          "fun_configs": ["PolytopeCSC_intersectCoord=cached_b_ref",
                          "PolytopeCSC_intersectCoord=cached_b_vec",
                          #"PolytopeCSC_intersectCoord=cached_b_vec_vec_nan", #same as nogather
                          "PolytopeCSC_intersectCoord=cached_b_vec_nogather",
                          #"PolytopeCSC_intersectCoord=cached_b_vec_nan", #same as nogather
                          "PolytopeCSC_intersectCoord=cached_b_vec_nan_inv",
                          #"PolytopeCSC_intersectCoord=cached_b_vec_inl", #same as vec
                          #"PolytopeCSC_intersectCoord=cached_b_vec_inl_2accs", #same as vec
          ],
          "run_configs": ["r=100000,polytopeType=2,intersect=intersectCoord_only"],
          "input_configs": [("generator", intersectSparseBodies)]
       },
    ],
    "xoption": ("generator", intersectSparseDims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
   },

   {"name": "csc_intersect_only_cross_cacheb",
    "executable": "benchmark_intersect",
    "config": [
       {
          "const_configs": [],
          "fun_configs": ["PolytopeCSC_intersectCoord=cached_b_ref",
                          "PolytopeCSC_intersectCoord=cached_b_vec",
                          "PolytopeCSC_intersectCoord=cached_b_vec_vec_nan",
                          "PolytopeCSC_intersectCoord=cached_b_vec_nogather",
                          "PolytopeCSC_intersectCoord=cached_b_vec_nan",
                          "PolytopeCSC_intersectCoord=cached_b_vec_nan_inv",
                          "PolytopeCSC_intersectCoord=cached_b_vec_inl",
                          "PolytopeCSC_intersectCoord=cached_b_vec_inl_2accs",
          ],
          "run_configs": ["r=100000,polytopeType=2,intersect=intersectCoord_only"],
          "input_configs": [("generator", crossBodies)]
       },
    ],
    "xoption": ("generator", crossDims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
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
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
   },

   {"name": "csc_readtest_dense",
    "executable": "benchmark_intersect",
    "config": [
       {
          "const_configs": [],
          "fun_configs": ["PolytopeCSC_intersectCoord=cached_vec_onlyread",
                          "PolytopeCSC_intersectCoord=cached_b_vec_nan_inv",
          ],
          "run_configs": ["r=100000,polytopeType=2,intersect=intersectCoord_only"],
          "input_configs": [("generator", intersectbodies)]
       },
    ],
    "xoption": ("generator", intersectdims),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["dim", "dim", "dim", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
   },

   ###################################################################################################
   #CSC END
   ###################################################################################################

   {"name": "get_rd_0_1",
    "executable": "benchmark_randomness",
    "config": [       
       {
          "const_configs": [],
          "fun_configs": rd_0_1,
          "run_configs": ["r=1000,rand_val_t=random_double_0_1"],
          "input_configs": [("i", [2**i for i in range(10,16)])],
       }
    ],
    "xoption": ("d", {str(2**i): str(2**i) for i in range(10,16)}),
    "title": ["Runtime Comparison", "Performance comparison", "I/O comparison", "Roofline measurements"],
    "xlabel": ["#doubles", "#doubles", "#doubles", "Operational Intensity [Flops/Byte]"],
    "ylabel": ["cycles(mean)", "flops/cylce(mean)", "bytes/cylce(mean)", "Performance [Flops/Cycle]"],
    "perf_roofs": [],
    "mem_roofs": []
   },
]
