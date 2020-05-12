#include <iostream>

extern "C" { // must be included C stlye
#include "../../src/volume/volume.h"
}

#include "../../src/volume/volume_helper.hpp"

#include "../../src/util/cli.hpp"
#include "../../src/util/cli_functions.hpp"

void test_box_inside(const int n, Body_T* type, void* body) {
   FT* x = (FT*)(aligned_alloc(32, n*sizeof(FT)));
   FT* d = (FT*)(aligned_alloc(32, n*sizeof(FT)));
   {
      for(int i=0;i<n;i++) {x[i]=0;}
      assert(type->inside(body, x));
   }
   {
      for(int i=0;i<n;i++) {x[i]=(i==2)*3;}
      assert(!type->inside(body, x));
   }
   {
      for(int i=0;i<n;i++) {x[i]=2;}
      assert(type->inside(body, x));
   }
   {
      for(int i=0;i<n;i++) {x[i]=1-2*(i%2==0);}
      assert(type->inside(body, x));
   }
   {
      for(int i=0;i<n;i++) {x[i]=-(i==0)*3.0;}
      assert(!type->inside(body, x));
   }
   free(x);
   free(d);
}

void test_box_intersect(const int n, Body_T* type, void* box) {
   FT* x = (FT*)(aligned_alloc(32, n*sizeof(FT)));
   FT* d = (FT*)(aligned_alloc(32, n*sizeof(FT)));

   // Test Polytope_T.intersect:
   {
      for(int i=0;i<n;i++) {x[i]=0; d[i]=-0.1 - (i==0)*0.9 + (i==1)*0.8;}
      FT t0,t1;
      type->intersect(box, x, d, &t0, &t1);
      assert(t0==-2.0 && t1==2.0);
   }
   {
      for(int i=0;i<n;i++) {x[i]=(i==2); d[i]=(i==2)*-1;}
      FT t0,t1;
      type->intersect(box, x, d, &t0, &t1);
      assert(t0==-1.0 && t1==3.0);
   }
   {
      for(int i=0;i<n;i++) {x[i]=(i==2)+(i==3)*-1.5; d[i]=0.1 + 0.4*(i==3);}
      FT t0,t1;
      type->intersect(box, x, d, &t0, &t1);
      assert(t0==-1.0 && t1==7.0);
   }
   {
      for(int i=0;i<n;i++) {x[i]=(i<2)*1.5; d[i]=1*(i==0) + -1*(i==1);}
      FT t0,t1;
      type->intersect(box, x, d, &t0, &t1);
      assert(t0==-0.5 && t1==0.5);
   }
   free(x);
   free(d);
}

void test_box_intersectCoord(const int n, Body_T* type, void* box) {
   FT* x = (FT*)(aligned_alloc(32, n*sizeof(FT)));
   FT* d = (FT*)(aligned_alloc(32, n*sizeof(FT)));

   // Test Polytope_T.intersectCoord
   void* cache = aligned_alloc(32, type->cacheAlloc(box));
   {
      for(int i=0;i<n;i++) {x[i]=0;}
      type->cacheReset(box,x,cache);
      for(int d=0;d<4;d++) {
         FT t0,t1;
         type->intersectCoord(box, x, 0, &t0, &t1, cache);
         assert(t0==-2.0 && t1==2.0);
      }
   }
   {
      for(int i=0;i<n;i++) {x[i]=(i==2);}
      type->cacheReset(box,x,cache);
      FT t0,t1;
      type->intersectCoord(box, x, 2, &t0, &t1, cache);
      assert(t0==-3.0 && t1==1.0);
      type->intersectCoord(box, x, 0, &t0, &t1, cache);
      assert(t0==-2.0 && t1==2.0);
      type->intersectCoord(box, x, 1, &t0, &t1, cache);
      assert(t0==-2.0 && t1==2.0);
      type->intersectCoord(box, x, 3, &t0, &t1, cache);
      assert(t0==-2.0 && t1==2.0);
   }
   free(cache);
   free(d);
   free(x);
}

void test_ellipsoid_intersectCoord(const int n, Body_T* type, Ellipsoid* e) {
   FT* x = (FT*)(aligned_alloc(32, n*sizeof(FT)));
   int d = prng_get_random_int_in_range(0,n-1);
   for(int i=0; i<n; i++) {
      x[i] = e->a[i] + prng_get_random_double_in_range(-1.0/n,1.0/n);
   }
      
   void* cache = aligned_alloc(32, Ellipsoid_T.cacheAlloc(e));
   Ellipsoid_T.cacheReset(e,x,cache);
   FT t0,t1;
   Ellipsoid_T.intersectCoord(e, x, d, &t0, &t1, cache);
      
   assert(t0 < t1);
      
   // test if points are really on ellipse:
   FT x0[n];
   FT x1[n];
   for(int i=0; i<n; i++) {
      x0[i] = x[i] + t0*(i==d);
      x1[i] = x[i] + t1*(i==d);
   }
   FT eval0 = Ellipsoid_eval(e, x0);
   FT eval1 = Ellipsoid_eval(e, x1);
   assert(std::abs(eval0-1) < 0.000001);
   assert(std::abs(eval1-1) < 0.000001);

   free(cache);
   free(x);
}


void test_box_cutOracle(const int n, Body_T* type, void* box) {
   FT* v = (FT*)(aligned_alloc(32, n*sizeof(FT)));
   FT c;

   Ellipsoid* e = Ellipsoid_new_with_T(n); // simple sphere
   for(int i=0; i<n; i++) {
      e->a[i] = prng_get_random_double_in_range(-0.1,0.1);
      FT* Ai = Ellipsoid_get_Ai(e,i);
      FT* Ti = Ellipsoid_get_Ai(e,i);
      FT r = prng_get_random_double_in_range(1.9*n,2.2*n);
      Ai[i] = 1.0 / (r*r);
      Ti[i] = (r*r);
   }
   
   {// fully inside inner ellipsoid:
      bool doCut = type->shallowCutOracle(box, e, v, &c);
      assert(!doCut && "center of ellipsoid");
   }

   for(int i=0;i<n;i++) {// center outside polytope
      for(int j=0;j<n;j++) { e->a[j] = (i==j)*(-3.0)*n + prng_get_random_double_in_range(-0.5,0.5); }

      bool doCut = type->shallowCutOracle(box, e, v, &c);
      assert(doCut && "outer ellipsoid");
      assert(c==1.0);
      for(int j=0;j<n;j++) { assert(v[j]==(i==j)*-1.0);}
   }

   for(int i=0;i<n;i++) {// center inside, but violate inner ellipsoid
      for(int j=0;j<n;j++) { e->a[j] = (i==j)*(-3.0) + prng_get_random_double_in_range(-0.1,0.1); }

      bool doCut = type->shallowCutOracle(box, e, v, &c);
      assert(doCut && "inner ellipsoid");
      assert(c==1.0);
      for(int j=0;j<n;j++) { assert(v[j]==(i==j)*-1.0);}
   }
   
   Ellipsoid_free(e);
   free(v);
}

int main(int argc, char** argv) {
   CLI cli(argc,argv,"test_volume_basics");
   CLIFunctionsVolume cliFun(cli);
  
   cliFun.preParse();
   if (!cli.parse()) {return -1;}
   cliFun.postParse();
   
   // -------------------------------- start tests
   
   // -------- dotProduct:
   {
      auto o = dynamic_cast<CLIF_Option<dotProduct_f_t>*>(cliFun.getOption("dotProduct"));
      for(auto it : o->fmap) {
         std::cout << "Test dotProduct " << it.first << " - " << it.second.second << std::endl;
         dotProduct = it.second.first;
         
         FT* u = (FT*)(aligned_alloc(32, 20*sizeof(FT)));
         FT* v = (FT*)(aligned_alloc(32, 20*sizeof(FT)));
         
         for(int i=0;i<20;i++) {
            u[i] = prng_get_random_double_in_range(1.1,2.0);
            v[i] = prng_get_random_double_in_range(1.1,2.0);
         }
         for(int n=1;n<20;n++){
            FT dot = dotProduct(u,v,n);
            
	    FT dotRef = 0;
	    for(int i=0;i<n;i++) {dotRef += u[i]*v[i];}
	    assert(std::abs(dot-dotRef) < 0.000001);
	 }
         free(u);
         free(v);
      }
   }
   // -------- squaredNorm:
   {
      auto o = dynamic_cast<CLIF_Option<squaredNorm_f_t>*>(cliFun.getOption("squaredNorm"));
      for(auto it : o->fmap) {
         std::cout << "Test squaredNorm " << it.first << " - " << it.second.second << std::endl;
         squaredNorm = it.second.first;
         
         FT* v = (FT*)(aligned_alloc(32, 20*sizeof(FT)));
         
         for(int i=0;i<20;i++) {
            v[i] = prng_get_random_double_in_range(1.1,2.0);
         }
         for(int n=1;n<20;n++){
            FT dot = squaredNorm(v,n);
            
	    FT dotRef = 0;
	    for(int i=0;i<n;i++) {dotRef += v[i]*v[i];}
	    assert(std::abs(dot-dotRef) < 0.000001);
	 }
         free(v);
      }
   }
 
   // --------------------------------- Polytope:
   auto o = dynamic_cast<CLIF_Option<intersectCoord_f_t>*>(cliFun.getOption("Polytope_intersectCoord"));
   for(auto it : o->fmap) {
      Polytope_T.intersectCoord = it.second.first;
      std::cout << "Test Polytope for intersectCoord " << it.first << " - " << it.second.second << std::endl;

      // Generate new polytope box, n dim, 2 radius
      for(int n=4;n<20;n++) {
	 std::cout << "test for n="<<n<<"\n";
         Polytope* box = Polytope_new_box(n,2);

         test_box_inside(n, &Polytope_T, box);
         test_box_intersect(n, &Polytope_T, box);
         test_box_intersectCoord(n, &Polytope_T, box);
         
         Polytope_free(box);
      }
   }

   auto oT = dynamic_cast<CLIF_Option<intersectCoord_f_t>*>(cliFun.getOption("PolytopeT_intersectCoord"));
   for(auto it : oT->fmap) {
      PolytopeT_T.intersectCoord = it.second.first;
      std::cout << "Test PolytopeT for intersectCoord " << it.first << " - " << it.second.second << std::endl;

      // Generate new polytope box, n dim, 2 radius
      for(int n=4;n<20;n++) {
	 std::cout << "test for n="<<n<<"\n";
         PolytopeT* box = PolytopeT_new_box(n,2);

         test_box_inside(n, &PolytopeT_T, box);
         test_box_intersect(n, &PolytopeT_T, box);
         test_box_intersectCoord(n, &PolytopeT_T, box);
         
         PolytopeT_free(box);
      }
   }


   auto oCSC = dynamic_cast<CLIF_Option<intersectCoord_f_t>*>(cliFun.getOption("PolytopeCSC_intersectCoord"));
   for (auto it : oCSC->fmap){
       // test PolytopeCSC
       std::cout << "Test PolytopeCSC for intersectCoord " << it.first << " - " << it.second.second << std::endl;

       // Generate new polytope box, n dim, 2 radius
       for(int n=4;n<20;n++) {
	  std::cout << "test for n="<<n<<"\n";
          Polytope* boxx = Polytope_new_box(n,2);
          PolytopeCSC *box = Polytope_to_PolytopeCSC(boxx);

          test_box_inside(n, &PolytopeCSC_T, box);
          test_box_intersect(n, &PolytopeCSC_T, box);
          test_box_intersectCoord(n, &PolytopeCSC_T, box);

          Polytope_free(boxx);
          PolytopeCSC_free(box);
       }
   }
   
   /*
   auto oJIT = dynamic_cast<CLIF_Option<intersectCoord_f_t>*>(cliFun.getOption("PolytopeJIT_intersectCoord"));
   for (auto it : oJIT->fmap){
       // test PolytopeJIT
       std::cout << "Test PolytopeJIT for intersectCoord " << it.first << " - " << it.second.second << std::endl;

       // Generate new polytope box, n dim, 2 radius
       for(int n=4;n<20;n++) {
	  std::cout << "test for n="<<n<<"\n";
          Polytope* boxx = Polytope_new_box(n,2);
          PolytopeJIT *box = Polytope_to_PolytopeJIT(boxx);

          test_box_inside(n, &PolytopeJIT_T, box);
          test_box_intersect(n, &PolytopeJIT_T, box);
          test_box_intersectCoord(n, &PolytopeJIT_T, box);

          Polytope_free(boxx);
          PolytopeJIT_free(box);
      }
   }
   */
   

   // Check ball volume:
   assert(std::abs(Ball_volume(3,1.0) - 4.189) <= 0.01);
   assert(std::abs(Ball_volume(4,1.0) - 4.935) <= 0.01);
   assert(std::abs(Ball_volume(10,1.0) - 2.550) <= 0.01);
   assert(std::abs(Ball_volume(11,1.0) - 1.884) <= 0.01);

   // ----------------- Ellipsoid
   { // test Ellipsoid_eval and Ellipsoid_normal Ellipsoid_T.inside
      const int n = 3;
      Ellipsoid* e = Ellipsoid_new(n); // simple sphere
      FT* normal = (FT*)(aligned_alloc(32, n*sizeof(FT)));
      FT* x = (FT*)(aligned_alloc(32, n*sizeof(FT)));
      {
         for(int i=0;i<n;i++) {x[i]=0;}
	 FT eval = Ellipsoid_eval(e, x);
	 assert(eval == 0);
	 assert(Ellipsoid_T.inside(e,x));
         Ellipsoid_normal(e, x, normal);
	 assert(normal[0]==0 && normal[1]==0 && normal[2]==0);
      }
      {
         for(int i=0;i<n;i++) {x[i]=(i==2);}
	 FT eval = Ellipsoid_eval(e, x);
	 assert(eval == 1.0);
	 assert(Ellipsoid_T.inside(e,x));
         Ellipsoid_normal(e, x, normal);
	 assert(normal[0]==0 && normal[1]==0 && normal[2]==2);
      }
      {
         for(int i=0;i<n;i++) {x[i]=(i==1)*2;}
	 FT eval = Ellipsoid_eval(e, x);
	 assert(eval == 4.0);
	 assert(!Ellipsoid_T.inside(e,x));
         Ellipsoid_normal(e, x, normal);
	 assert(normal[0]==0 && normal[1]==4 && normal[2]==0);
      } 
      {
         for(int i=0;i<n;i++) {x[i]=(i==0)*3 + (i==1)*4;}
	 FT eval = Ellipsoid_eval(e, x);
	 assert(eval == 25.0);
	 assert(!Ellipsoid_T.inside(e,x));
         Ellipsoid_normal(e, x, normal);
	 assert(normal[0]==6 && normal[1]==8 && normal[2]==0);
      }
      free(x);
      free(normal);
      Ellipsoid_free(e);
   }
   { // test Ellipsoid_T.inside
      const int n = 20;
      Ellipsoid* e = Ellipsoid_new(n); // simple sphere
      for(int i=0; i<n; i++) {
         FT* Ai = Ellipsoid_get_Ai(e,i);
         Ai[i] = (i+1)*(i+1);
      }
      FT* x = (FT*)(aligned_alloc(32, n*sizeof(FT)));
      {
         for(int i=0; i<n; i++) {x[i] = 0;}
	 assert(Ellipsoid_T.inside(e,x));
      }
      for(int j=0;j<n;j++){
         for(int i=0; i<n; i++) {x[i] = (i==j)*(1.0/n);}
	 assert(Ellipsoid_T.inside(e,x));
         for(int i=0; i<n; i++) {x[i] = (i==j)*(1.0/(i+1.0));}
	 assert(Ellipsoid_T.inside(e,x));
         for(int i=0; i<n; i++) {x[i] = (i==j)*(1.1/(i+1.0));}
	 assert(!Ellipsoid_T.inside(e,x));
      }
      free(x);
      Ellipsoid_free(e);
   }
   for(int t=0; t<100; t++) {// test Ellipsoid_T.intersect
      const int n = 20;
      Ellipsoid* e = Ellipsoid_new(n); // simple sphere
      for(int i=0; i<n; i++) {
         e->a[i] = prng_get_random_double_in_range(-10,10);
         FT* Ai = Ellipsoid_get_Ai(e,i);
         Ai[i] = prng_get_random_double_in_range(1.1,2.0);
      }
      
      FT* x = (FT*)(aligned_alloc(32, n*sizeof(FT)));
      FT* d = (FT*)(aligned_alloc(32, n*sizeof(FT)));
      for(int i=0; i<n; i++) {
         x[i] = e->a[i] + prng_get_random_double_in_range(-1.0/n,1.0/n);
         d[i] = prng_get_random_double_normal();
      }
      
      FT t0,t1;
      Ellipsoid_T.intersect(e, x, d, &t0, &t1);
      
      assert(t0 < t1 && std::abs(t1-t0) > 0.1);
      
      // test if points are really on ellipse:
      FT x0[n];
      FT x1[n];
      for(int i=0; i<n; i++) {
         x0[i] = x[i] + t0*d[i];
         x1[i] = x[i] + t1*d[i];
      }
      FT eval0 = Ellipsoid_eval(e, x0);
      FT eval1 = Ellipsoid_eval(e, x1);
      assert(std::abs(eval0-1) < 0.000001);
      assert(std::abs(eval1-1) < 0.000001);

      free(d);
      free(x);
      Ellipsoid_T.free(e);
   }

   auto oE = dynamic_cast<CLIF_Option<intersectCoord_f_t>*>(cliFun.getOption("Ellipsoid_intersectCoord"));
   for(auto it : oE->fmap) {
      Ellipsoid_T.intersectCoord = it.second.first;
      std::cout << "Test Ellipsoid for intersectCoord " << it.first << " - " << it.second.second << std::endl;

      for (int t = 0; t < 100; t++) {
         // Generate new ellipsoid box, n dim
         const int n = 20;
         Ellipsoid* e = Ellipsoid_new(n); // simple sphere
         for(int i=0; i<n; i++) {
            e->a[i] = prng_get_random_double_in_range(-10,10);
            FT* Ai = Ellipsoid_get_Ai(e,i);
            Ai[i] = prng_get_random_double_in_range(1.1,2.0);
         }

         test_ellipsoid_intersectCoord(n, &Ellipsoid_T, e);
         
         Ellipsoid_T.free(e);
      }
   }

   // Test Ellipsoid cache update
   auto oEc = dynamic_cast<CLIF_Option<cacheUpdateCoord_f_t>*>(cliFun.getOption("Ellipsoid_cacheUpdateCoord"));
   for(auto it : oEc->fmap) {
      Ellipsoid_T.cacheUpdateCoord = it.second.first;
      std::cout << "Test Ellipsoid for cacheUpdateCoord " << it.first << " - " << it.second.second << std::endl;

      for (int t = 0; t < 100; t++) {
         // Generate new ellipsoid box, n dim
         const int n = 20;
         Ellipsoid* e = Ellipsoid_new(n); // simple sphere
         for(int i=0; i<n; i++) {
            e->a[i] = prng_get_random_double_in_range(-10,10);
            FT* Ai = Ellipsoid_get_Ai(e,i);
            Ai[i] = prng_get_random_double_in_range(1.1,2.0);
         }

         test_ellipsoid_intersectCoord(n, &Ellipsoid_T, e);
         
         Ellipsoid_T.free(e);
      }
   }

   { // test Ellipsoid_project
      Ellipsoid* e = Ellipsoid_new(3); // simple sphere
      e->a[0] = 5.0;
      e->a[1] = -10.0;
      e->a[2] = 0.5;
      
      FT* A0 = Ellipsoid_get_Ai(e,0);
      FT* A1 = Ellipsoid_get_Ai(e,1);
      FT* A2 = Ellipsoid_get_Ai(e,2);
      A0[0] = 2.0;
      A1[1] = 1.0;
      A2[2] = 0.5;

      FT* x = (FT*)(aligned_alloc(32, 3*sizeof(FT)));
      for(int t=0; t<100; t++) {
         x[0] = std::fmod(t,13);
	 x[1] = std::fmod(t,11) - 12;
	 x[2] = std::fmod(t,0.7);
         FT r = prng_get_random_double_in_range(0.1,20);
	 Ellipsoid_project(e, r, x);
         FT eval = Ellipsoid_eval(e, x);
         assert(std::abs(eval-r) < 0.00001);
      }
      free(x);
   }
   
   for(int t=0; t<100; t++) {// test Ellipsoid_minimize
      const int n = 20;
      Ellipsoid* e1 = Ellipsoid_new(n); // simple sphere
      Ellipsoid* e2 = Ellipsoid_new(n); // simple sphere
      for(int i=0; i<n; i++) {
         e1->a[i] = prng_get_random_double_in_range(-10,10);
         e2->a[i] = prng_get_random_double_in_range(-10,10);
         FT* A1i = Ellipsoid_get_Ai(e1,i);
         A1i[i] = prng_get_random_double_in_range(0.1,20);
         FT* A2i = Ellipsoid_get_Ai(e2,i);
         A2i[i] = prng_get_random_double_in_range(0.1,20);
      }
      
      {
         FT* x = (FT*)(aligned_alloc(32, n*sizeof(FT)));
         FT* n1 = (FT*)(aligned_alloc(32, n*sizeof(FT)));
         FT* n2 = (FT*)(aligned_alloc(32, n*sizeof(FT)));
         FT* step = (FT*)(aligned_alloc(32, n*sizeof(FT)));
         for(int i=0; i<n; i++) {
	    x[i] = e2->a[i] + prng_get_random_double_normal();
	 }

	 FT r = prng_get_random_double_in_range(0.1,20);// change size of e2
	 Ellipsoid_minimize(e2,r,e1,x); // test
         FT eval = Ellipsoid_eval(e2,x);
	 std::cout << r << " " << eval << "\n";
	 assert(std::abs(eval-r) < 0.000001);
	 Ellipsoid_normal(e1, x, n1);
         Ellipsoid_normal(e2, x, n2);
	 FT n1_2 = dotProduct(n1,n1, n);
	 FT dot = dotProduct(n1,n2, n);
	 
         for(int i=0;i<n;i++) { step[i] = n2[i] - n1[i] * dot / n1_2;}
	 FT step2 = dotProduct(step,step,n);
	 assert(std::abs(step2) < 0.0001);
	 
	 free(step);
	 free(n2);
	 free(n1);
	 free(x);
      }
      

      Ellipsoid_free(e1);
      Ellipsoid_free(e2);
   }
   
   // --------------------------------------------------- Shallow Cut Oracles
   {// Polytope_T.shallowCutOracle
      const int n = 20;
      Polytope* box = Polytope_new_box(n,1.0);
      
      test_box_cutOracle(n, &Polytope_T, box);
      
      Polytope_free(box);
   }
   {// PolytopeT_T.shallowCutOracle
      const int n = 20;
      PolytopeT* box = PolytopeT_new_box(n,1.0);
      
      test_box_cutOracle(n, &PolytopeT_T, box);
      
      PolytopeT_free(box);
   }
   { // PolytopeCSC_T.shallowCutOracle
      const int n = 20;
      Polytope* bbox = Polytope_new_box(n,1.0);
      PolytopeCSC *box = Polytope_to_PolytopeCSC(bbox);
      
      test_box_cutOracle(n, &PolytopeCSC_T, box);
      
      PolytopeT_free(bbox);
      PolytopeCSC_free(box);
   }
   {// Ellipsoid_T.shallowCutOracle
      std::cout << "Ellipsoid oracle:\n";
      const int n = 20;
      Ellipsoid* body = Ellipsoid_new(n); // simple sphere
      
      FT* v = (FT*)(aligned_alloc(32, n*sizeof(FT)));
      FT c;
      

      Ellipsoid* e = Ellipsoid_new_with_T(n); // simple sphere
      for(int i=0; i<n; i++) {
         e->a[i] = prng_get_random_double_in_range(-0.001,0.001);
         FT* Ai = Ellipsoid_get_Ai(e,i);
         FT* Ti = Ellipsoid_get_Ti(e,i);
         FT r = prng_get_random_double_in_range(1.9*n,1.9*n); // make inner ellipsoid smaller than sphere!
         Ai[i] = 1.0 / (r*r);
         Ti[i] = (r*r);
      }
      
      std::cout << "fully inside:\n";
      {// fully inside inner ellipsoid:
         // sanity check:
         assert(Ellipsoid_T.inside(body, e->a));
         
         bool doCut = Ellipsoid_T.shallowCutOracle(body, e, v, &c);
         assert(!doCut && "center of ellipsoid");
      }

      std::cout << "center outside:\n";
      for(int i=0;i<n;i++) {// center outside polytope
         for(int j=0;j<n;j++) { e->a[j] = (i==j)*(-3.0)*n + prng_get_random_double_in_range(-0.5,0.5); }

	 // sanity check:
	 assert(!Ellipsoid_T.inside(body, e->a));
        
         bool doCut = Ellipsoid_T.shallowCutOracle(body, e, v, &c);
         assert(doCut && "outer ellipsoid");
         
         // check that center of body is in halfspace:
         FT dot = dotProduct(v, body->a, n);
         assert(dot <= c);
      }
      
      std::cout << "center sandwich:\n";
      for(int i=0;i<n;i++) {// center inside, but violate inner ellipsoid
         for(int j=0;j<n;j++) { e->a[j] = (i==j)*(-0.3) + prng_get_random_double_in_range(-0.01,0.01); }
 
	 // sanity check:
	 assert(Ellipsoid_T.inside(body, e->a));
         
	 bool doCut = Ellipsoid_T.shallowCutOracle(body, e, v, &c);
         assert(doCut && "inner ellipsoid");
      	 
         // check that center of body is in halfspace:
         FT dot = dotProduct(v, body->a, n);
         assert(dot <= c);
      }

      Ellipsoid_free(e);
      free(v);
      Ellipsoid_free(body);
   }

   // ---------------------------------------- Matrix

   for(int t=0;t<100;t++){// Matrix_L_solve:
      const int n = 20;
      Matrix* L = Matrix_new(n,n);
      FT* x = (FT*)(aligned_alloc(32, n*sizeof(FT)));
      FT* b = (FT*)(aligned_alloc(32, n*sizeof(FT)));
      // init L matrix:
      for(int i=0;i<n;i++) {
         FT* Li = Matrix_get_row(L,i);
	 Li[i] = prng_get_random_double_in_range(5,10);// strong diagonal
	 for(int j=0;j<i;j++) {
	    Li[j] = prng_get_random_double_in_range(0.1,0.5);
	 }
      }
      
      for(int i=0;i<n;i++) {
         for(int j=0;j<n;j++) {b[j]=(i==j);} // for each unit vector
	 Matrix_L_solve(L,x,b);

	 for(int j=0;j<n;j++) {// check by multiplying
            FT* Li = Matrix_get_row(L,j);
	    FT dot = dotProduct(Li,x, n);
	    assert(std::abs(dot - b[j])<0.000001);
	 }
      } 

      Matrix_free(L);
      free(x);
      free(b);
   }


   //Polytope *P, *Q;
   //FT vol1, vol2;
   //read_vinci("../polyvest/examples/rh_10_20.ine", &P, &vol1);
   //Polytope_print((void *) P);
   //cout << vol1 << "\n";

   //read_vinci("../polyvest/examples/cc_8_10.ine", &Q, &vol2);
   //Polytope_print((void *) Q);
   //cout << vol2 << "\n";


   
   // -------------------------------- end tests

   #ifdef NDEBUG
   std::cout<< "WARNING: DEBUG DISABLED!\n";
   #else
   std::cout<< "TESTS COMPLETE.\n";
   #endif
}

