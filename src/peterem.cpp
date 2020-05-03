// Please ignore this file, it is only used for small WIP tests

#include "peterem.hpp"


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


   cliFun.preParse();
   if (!cli.parse()) {return -1;}
   cliFun.postParse();
   

   Solved_Body* solved_body = solved_body_generator()->get(generator,polytopeTranspose);
   solved_body->print();
   const int n = solved_body->n;

   FT* p = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
   
   Ellipsoid* e = Ellipsoid_new(n);

   auto f = [&](FT x, FT y, FT &r, FT &g, FT &b) {
      FT xx = x*4*n-2*n;
      FT yy = y*4*n-2*n;
      for(int i=0;i<n;i++) {p[i] = (i==0)*xx + (i==1)*yy;}
     
      FT eval = Ellipsoid_eval(e,p);
      r = 0;
      g = eval*(eval<=1.0);
      b = eval/(4*n*n)*(eval <= 4*n*n);
      

      auto &type = solved_body->type;
      auto &body = solved_body->body;
      auto &bcount = solved_body->bcount;
      for(int c=0;c<bcount;c++) {
         if(type[c]->inside(body[c],p)) {
	    r += 0.5;
	 }
      }
   };

   {// ---------------------- IMG
      int width = 800;
      int height = 800;
      EVP::Image_BMP img(height,width);
      for(int i=0; i<height; i++){
          for(int j=0; j<width; j++){
	      FT rr,gg,bb;
	      f((double)i/height, (double)j/width, rr,gg,bb);
              int r = (unsigned char)(rr*255); ///red
              int g = (unsigned char)(gg*255); ///green
              int b = (unsigned char)(bb*255); ///blue
              img.set(j,i,r,g,b);
	  }
      }
      img.toFile("peterem_out.bmp");
   }// ---------------------- END IMG

   #ifdef NDEBUG
   std::cout<< "## WARNING: DEBUG DISABLED!\n";
   #else
   std::cout<< "## TESTS COMPLETE.\n";
   #endif
}




