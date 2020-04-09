#!/usr/bin/env python3

import sys
import os
import subprocess
import pprint
from threading import Timer

# --------------------------------- ADD YOUR TEST HERE
TESTS = ["volume/test_volume_basics",
         "volume/test_volume_estimate",
         "preprocess/test_preprocess",
        ];


# ---------------------------- Parse python args
DO_TESTS = TESTS;
if(len(sys.argv)>1):
   DO_TESTS = [];
   
   for i in range(1,len(sys.argv)):
      test = sys.argv[i]
      if(test in TESTS):
         DO_TESTS.append(test);
      else:
         print("ERROR: test '{}' is not available!".format(test));
         print("  list of available tests:");
         for t in TESTS:
            print("     ",t);
         sys.exit(1);

# ------ iterate over tests chosen in DO_TESTS, subset of TESTS

def run_test(test):
   print("# Running Test '{}'...".format(test));
   myenv = os.environ;
   #myenv["OMP_NUM_THREADS"] = str(nproc); # change env
   proc = subprocess.Popen((sys.path[0]+"/"+test,), stdout=subprocess.PIPE, stderr=subprocess.PIPE, env = myenv);
   
   isTimeout = False;
   try:
      outs, errs = proc.communicate(timeout=20)
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
      print("# FAIL: '{}'".format(test));
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
