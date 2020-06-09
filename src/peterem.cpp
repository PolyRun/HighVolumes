// Please ignore this file, it is only used for small WIP tests

#include "peterem.hpp"

void optimize_test(const std::string& generator) {
   Solved_Body* solved_body = solved_body_generator()->get(generator,false);
   solved_body->print();
   const int n = solved_body->n;

   if(solved_body->bcount ==1 && solved_body->type[0] == &Polytope_T) {
      Polytope* p = (Polytope*)solved_body->body[0];
   
      Polytope* q = optimize_polytope(p);
      std::cout << "\n## Optimized: \n\n";
      Polytope_T.print(q);
   } else {
      std::cout << "\n could not optimize, not a single polytope!\n";
   }
}


int main(int argc, char** argv) {
   CLI cli(argc,argv,"peterem");
   CLIFunctionsVolume cliFun(cli);
   
   cliFun.claimOpt('b',"body generation configuration");

   std::string generator = "cube";
   auto &gen_map = solved_body_generator()->gen_map();
   cliFun.add(new CLIF_Option<std::string>(&generator,'b',"generator","cube_r1.0_10", gen_map));
   
   bool polytopeTranspose = false;
   cliFun.add(new CLIF_Option<bool>(&polytopeTranspose,'b',"polytopeTranspose","false", {
                                                    {"false",{false, "Polytope format / rows"}},
       					     {"true",{true, "PolytopeT format / columns"}} }));

   
   int v1;
   double v2;
   cliFun.add(new CLIF_DoubleOption<int,double>(&v1,&v2,'e',"choice","c1", {
			   {"c1", {{1,1.1}, "desc c1"}},
			   {"c2", {{2,2.1}, "desc c2"}},
			   }));

   int v3;
   double v4;
   size_t v5;
   cliFun.add(new CLIF_TrippleOption<int,double,size_t>(&v3,&v4,&v5,'e',"choice2","c1", {
			   {"c1", {{1,{1.1, (size_t)0x1 << 35}}, "desc c1"}},
			   {"c2", {{1,{2.1, (size_t)0x2 << 35}}, "desc c2"}},
			   }));




   cliFun.preParse();
   if (!cli.parse()) {return -1;}
   cliFun.postParse();
   
   std::cout << "choice "<< v1 << " " << v2 << "\n";
   std::cout << "choice "<< v3 << " " << v4 << " " << v5 << "\n";
   
   auto a = ArbitraryExpNum_new(10);
   ArbitraryExpNum_print(a); printf(" - start\n");

   for(int i=0;i<400; i++) {
      a = ArbitraryExpNum_mul(a,1.0/11);
      ArbitraryExpNum_print(a); printf(" - %d\n",i);
   }
   for(int i=0;i<400; i++) {
      a = ArbitraryExpNum_mul(a,11);
      ArbitraryExpNum_print(a); printf(" - %d\n",i);
   }

   optimize_test(generator);

   Solved_Body* solved_body = solved_body_generator()->get(generator,polytopeTranspose);
   solved_body->print();
   const int n = solved_body->n;

   FT* p = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
   Ellipsoid* e = Ellipsoid_new(n);

   auto f = [&](FT x, FT y, FT z, FT &r, FT &g, FT &b) {
      FT xx = x*4*n-2*n;
      FT yy = y*4*n-2*n;
      FT zz = (z*4*n-2*n)/4.0;
      for(int i=0;i<n;i++) {p[i] = (i==0)*xx + (i==1)*yy + (i==2)*zz;}
     
      FT eval = Ellipsoid_eval(e,p);
      r = 0;
      g = eval*(eval<=1.0);
      b = eval/(4*n*n)*(eval <= 4*n*n);
      

      auto &type = solved_body->type;
      auto &body = solved_body->body;
      auto &bcount = solved_body->bcount;
      for(int c=0;c<bcount;c++) {
         if(type[c]->inside(body[c],p)) {
	    r+= 1.0/bcount;
	 }
      }
   };

   FT* dx = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
   FT* dx2 = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
   FT* dy = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
   FT* dy2 = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
   FT* dz = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
   FT* dz2 = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
   for(int i=0;i<n;i++) {dx[i] = (i==0);dy[i]=(i==1);dz[i]=(i==2);}
   FT* normal = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
   FT* p2 = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
   
   Matrix* R = Matrix_new(n,n); // rotation matrix
   for(int i=0;i<n;i++) {Matrix_set(R,i,i,1.0);}
   for(int i=0;i<n;i++) {
      for(int j=0;j<n;j++) {
         FT angle = prng_get_random_double_in_range(-0.2/n,0.2/n);//2*M_PI); 
         Matrix_rotate(R, i, j, angle);
      }
   }

   auto f2 = [&](FT x, FT y, FT z, FT &r, FT &g, FT &b) {
      FT xx = 4*x-2;//x*4*n-2*n;
      FT yy = 4*y-2;//y*4*n-2*n;
      FT zz = 0;
      for(int i=0;i<n;i++) {p[i] = dx[i]*xx + dy[i]*yy + dz[i]*zz;}
     
      FT eval = Ellipsoid_eval(e,p);
      r = 0;//2*n;
      g = 0;//eval*(eval<=1.0);
      b = 0;//eval/(4*n*n)*(eval <= 4*n*n);
      
      FT zzz = 2*n;

      auto &type = solved_body->type;
      auto &body = solved_body->body;
      auto &bcount = solved_body->bcount;
      for(int c=0;c<bcount;c++) {
	 FT t0,t1;
         type[c]->intersect(body[c], p, dz, &t0, &t1);
	 if(t0 < t1) {
	 //if(type[c]->inside(body[c],p)) {
	 //   FT t0,t1;
         //   type[c]->intersect(body[c], p, dz, &t0, &t1);
	    if(t1 < zzz) {
	       zzz = t1;
	       // update normal!
               for(int i=0;i<n;i++) {p2[i]=p[i] + dz[i]*(t1+0.001);}
               type[c]->normal(body[c], p2, normal);
	       FT n2 = squaredNorm(normal,n);
	       FT nInv = 1.0/sqrt(n2);
	       r = dotProduct(dx,normal,n)*nInv*0.5+0.5;
	       g = dotProduct(dy,normal,n)*nInv*0.5+0.5;
	       b = dotProduct(dz,normal,n)*nInv*0.5+0.5;
	    }
	 } else {
	    zzz = -1e20;
	 }
      }
      if(zzz < -1e19) {r = 0;g=0;b=0;}
   };

   {// ---------------------- IMG
      int width = 400;
      int height = 400;
      int depth = 200;
      EVP::Image_BMP img(height,width);
      for(int z=0;z<depth; z++) {
          Matrix_MVM(R,dx,dx2); std::swap(dx,dx2);
          Matrix_MVM(R,dy,dy2); std::swap(dy,dy2);
          Matrix_MVM(R,dz,dz2); std::swap(dz,dz2);
     	  for(int i=0; i<height; i++){
             for(int j=0; j<width; j++){
                 FT rr,gg,bb;
                 f2((double)i/height, (double)j/width, (double)z/depth, rr,gg,bb);
                 int r = (unsigned char)(rr*255); ///red
                 int g = (unsigned char)(gg*255); ///green
                 int b = (unsigned char)(bb*255); ///blue
                 img.set(j,i,r,g,b);
             }
         }
         std::string num = std::to_string(z);
	 while(num.size() < 5) {num = "0"+num;}
         img.toFile("out/peterem_out_"+num+".bmp");
      }
   }// ---------------------- END IMG

   #ifdef NDEBUG
   std::cout<< "## WARNING: DEBUG DISABLED!\n";
   #else
   std::cout<< "## TESTS COMPLETE.\n";
   #endif
}




