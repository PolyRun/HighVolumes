#!/usr/bin/env python3

import sys
import os
import subprocess
import pprint
from threading import Timer

# --------------------------------- ADD YOUR TEST HERE
TESTS = [
         {"name": "random/test_random",
          "timeout": 10},
         {"name": "jit/test_jit",
          "timeout": 10},
         {"name": "volume/test_volume_basics",
          "timeout": 50},
         {"name": "volume/test_volume_estimate",
          "timeout": 20},
         {"name": "volume/test_volume_examples",
          "timeout": 40},
         {"name": "preprocess/test_preprocess",
          "timeout": 30},
         {"name": "preprocess/test_init",
          "timeout": 10},
         {"name": "preprocess/unit_test",
          "timeout": 10},
         {"name": "matrix/matrix_unit_test",
          "timeout": 10},
        ];


# ---------------------------- Parse python args
DO_TESTS = TESTS;
if(len(sys.argv)>1):
   DO_TESTS = [];
   
   for i in range(1,len(sys.argv)):
      test = sys.argv[i]
      upd = list(filter(lambda t: t["name"] == test, TESTS))
      DO_TESTS += upd
      if len(upd) == 0:
         print("ERROR: test '{}' is not available!".format(test));
         print("  list of available tests:");
         for t in TESTS:
            print("     ",t);
         sys.exit(1);

# ------ iterate over tests chosen in DO_TESTS, subset of TESTS

def run_test(test):
   tname = test["name"]
   print("# Running Test '{}'...".format(tname));
   myenv = os.environ;
   #myenv["OMP_NUM_THREADS"] = str(nproc); # change env
   proc = subprocess.Popen((sys.path[0]+"/"+tname,), stdout=subprocess.PIPE, stderr=subprocess.PIPE, env = myenv);
   
   isTimeout = False;
   try:
      outs, errs = proc.communicate(timeout=test["timeout"])
   except subprocess.TimeoutExpired:
      proc.kill()
      outs, errs = proc.communicate()
      isTimeout = True;
   
   if(b"TESTS COMPLETE." in outs and errs == b"" and not isTimeout):
      print("# SUCCESS");
      return True;
   else:
      for o in outs.split(b"\n"):
         print("  > {}".format(o.decode("utf-8")));
      if(errs!=b""):
         for e in errs.split(b"\n"):
            print("  stderr: {}".format(e.decode("utf-8")));
      if(isTimeout):
         print("  TIMEOUT!");
      print("# FAIL: '{}'".format(tname));
      return False;

FAILED = []
for test in DO_TESTS:
   if(not run_test(test)):
      FAILED.append(test)
   print("");

n = len(DO_TESTS)
f = len(FAILED)
print("-"*50);
print("PASSED: {} of {}.".format(n-f,n))
for ff in FAILED:
   print(" FAILED: {}".format(ff))
print("-"*50);
