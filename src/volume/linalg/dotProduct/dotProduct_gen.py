#!/usr/bin/env python3
import math



# config:
cutOff_n = 100;
dotProduct = "auto2";

# ---------------------------------- emit
indent = 0;
def emit(s):
   global indent;
   print(" "*(3*indent) + s);


# --------------------------------- main
def emitUnvectorized():
   emit("FT sum = 0.0;");
   emit("for(int i=0; i<n; i++) {sum+= u[i]*v[i];}");
   emit("return sum;");

def emitVectorized(n):
   global indent;
  
   m = math.ceil(n / 4);
   rest = n % 4;
   emit("// m = "+str(m));
   emit("// rest = "+str(rest));

   # --- load:
   for i in range(m):
      if(i<m-1 or rest==0):
         emit("const __m256d u{} = _mm256_load_pd(u+{});".format(i,i*4));
         emit("const __m256d v{} = _mm256_load_pd(v+{});".format(i,i*4));
      else:
         on = "0x"+"F"*16
         emit("const __m256i mask = _mm256_set_epi64x({}, {}, {}, {});".format(
             0,
             on if (rest>=3) else 0,
             on if (rest>=2) else 0,
             on if (rest>=1) else 0));
         emit("const __m256d u{} = _mm256_maskload_pd(u+{},mask);".format(i,i*4));
         emit("const __m256d v{} = _mm256_maskload_pd(v+{},mask);".format(i,i*4));

   # --- compute:
   emit("__m256d sum = _mm256_mul_pd(u{},v{});".format(0,0));
   if(m==2):
      emit("sum = _mm256_fmadd_pd(u{},v{},sum);".format(1,1));
   elif(m>2):
      emit("__m256d sum2 = _mm256_mul_pd(u{},v{});".format(1,1));
      for i in range(2,m):
         x = "" if (i%2==0) else "2";
         emit("sum{} = _mm256_fmadd_pd(u{},v{},sum{});".format(x,i,i,x));
      emit("sum = _mm256_add_pd(sum,sum2);");
   
   # --- collapse:
   # can probably be improved
   emit("sum = _mm256_hadd_pd(sum,sum);");
   #emit("sum = _mm256_permute4x64_pd(sum, 0b00100010);");
   #emit("sum = _mm256_hadd_pd(sum,sum);");
 
   # --- return:
   emit("return sum[0]+sum[2];");

def emitScalar(n):
   global indent;

   terms = []
   for i in range(n):
      terms.append("u[{}]*v[{}]".format(i,i));
   terms = " + ".join(terms);
   emit("return {};".format(terms));

emit("// Code Generated by dotProduct_gen.py");
emit('#include "dotProduct.h"')
emit("FT dotProduct_"+dotProduct+"(const FT* u, const FT* v, const int n) {");
indent+=1
emit("switch(n) {");
indent+=1

emit("case 0: {"); # ------------------- 0
indent+=1

emit("return 0;");

emit("break;");
indent-=1
emit("}//case 0");

for i in range(1,cutOff_n+1):
   emit("case "+str(i)+": {");
   indent+=1
   
   if(True):
      emitScalar(i);
   else:
      emitVectorized(i);
   
   emit("break;");
   indent-=1
   emit("}//case "+str(i));


emit("default:{"); # ------------------- default
indent+=1

emitUnvectorized();

emit("break;");
indent-=1
emit("}//case default");



indent-=1
emit("}//switch");
indent-=1
emit("}//dotProduct_"+dotProduct);
