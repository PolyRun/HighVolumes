#!/usr/bin/env python3
import subprocess
import sys
import time

# switch polytopeType:
# 1 for our impl,
# 4 for PolyVest

for i in range(17,20):
   print("#-#-#-#-#-#-#-#-#-#-# Running {}".format(i));
   
   # benchmark/benchmark_app
   # -b "generator=birk_3,polytopeType=1,r=1"
   # -c "step_size=1600,verbose=1"
   # -f "rand_f=sr_rand,walk_f=walkCoord_8"

   start = time.time()
   proc = subprocess.Popen(
      [sys.path[0]+"/"+"/benchmark_app",
        '-b', 'generator=birk_{},polytopeType=1,r=1'.format(i),
        '-c', 'step_size=1600,verbose=1',
        '-f', "rand_f=sr_rand,walk_f=walkCoord_8",
        ],
      #stdout=subprocess.PIPE,
   );
   proc.wait()
   end = time.time()
   print("Time: {}".format(end - start))
   


