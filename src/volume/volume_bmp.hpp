// Header file to generate .BMP images from volume bodies


#ifndef HEADER_VOLUME_BMP_HPP
#define HEADER_VOLUME_BMP_HPP

#include "volume_helper.hpp"
#include "volume_examples.hpp"
#include "../util/image_bmp.hpp"

class Volume_BMP {
public:
   Volume_BMP(const std::string &path) {
      Solved_Body* solved_body = solved_body_generator()->get("2sphere_preprocessed_2",true);
      solved_body->print();
	
      
      int width = 400;
      int height = 300;
      EVP::Image_BMP img(height,width);
      for(int i=0; i<height; i++){
          for(int j=0; j<width; j++){
              int r = (unsigned char)((double)i/height*255); ///red
              int g = (unsigned char)((double)j/width*255); ///green
              int b = (unsigned char)(((double)i+j)/(height+width)*255); ///blue
              img.set(j,i,r,g,b);
	  }
      }
      img.toFile(path + "out.bmp");
   }
private:
};

#endif // HEADER_VOLUME_BMP_HPP


