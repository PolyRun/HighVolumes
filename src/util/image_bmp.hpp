// Helper file to generate .BMP image files

#ifndef HEADER_IMAGE_BMP_HPP
#define HEADER_IMAGE_BMP_HPP

#include <stdio.h>
#include <vector>

namespace EVP {
   const int bytesPerPixel = 3; /// red, green, blue
   const int fileHeaderSize = 14;
   const int infoHeaderSize = 40;

   class Image_BMP {
   public:
      Image_BMP(int height_, int width_) : height(height_),width(width_) {
         image.resize(height*width*bytesPerPixel,0);
      }

      void set(int w, int h, unsigned char r, unsigned char g, unsigned char b) {
         image[bytesPerPixel * (h*width + w) + 2] = r;
         image[bytesPerPixel * (h*width + w) + 1] = g;
         image[bytesPerPixel * (h*width + w) + 0] = b;
      }

      void toFile(const std::string &fname) {
         unsigned char padding[3] = {0, 0, 0};
         int paddingSize = (4 - (width*bytesPerPixel) % 4) % 4;
         
         auto fHeader = fileHeader(paddingSize);
         auto iHeader = infoHeader();
         
         FILE* imageFile = fopen(fname.c_str(), "wb");
         
         fwrite(fHeader.data(), 1, fileHeaderSize, imageFile);
         fwrite(iHeader.data(), 1, infoHeaderSize, imageFile);
         
         for(int i=0; i<height; i++){
             fwrite(image.data()+(i*width*bytesPerPixel), bytesPerPixel, width, imageFile);
             fwrite(padding, 1, paddingSize, imageFile);
         }
         
         fclose(imageFile);
      }

      std::vector<unsigned char> fileHeader(const int paddingSize) {
         int fileSize = fileHeaderSize + infoHeaderSize + (bytesPerPixel*width+paddingSize) * height;
         
	 std::vector<unsigned char> fileHeader = {
             0,0, /// signature
             0,0,0,0, /// image file size in bytes
             0,0,0,0, /// reserved
             0,0,0,0, /// start of pixel array
         };
         
         fileHeader[ 0] = (unsigned char)('B');
         fileHeader[ 1] = (unsigned char)('M');
         fileHeader[ 2] = (unsigned char)(fileSize    );
         fileHeader[ 3] = (unsigned char)(fileSize>> 8);
         fileHeader[ 4] = (unsigned char)(fileSize>>16);
         fileHeader[ 5] = (unsigned char)(fileSize>>24);
         fileHeader[10] = (unsigned char)(fileHeaderSize + infoHeaderSize);
         
         return fileHeader;
      }

      std::vector<unsigned char> infoHeader() {
	 std::vector<unsigned char> infoHeader = {
             0,0,0,0, /// header size
             0,0,0,0, /// image width
             0,0,0,0, /// image height
             0,0, /// number of color planes
             0,0, /// bits per pixel
             0,0,0,0, /// compression
             0,0,0,0, /// image size
             0,0,0,0, /// horizontal resolution
             0,0,0,0, /// vertical resolution
             0,0,0,0, /// colors in color table
             0,0,0,0, /// important color count
         };
         
         infoHeader[ 0] = (unsigned char)(infoHeaderSize);
         infoHeader[ 4] = (unsigned char)(width    );
         infoHeader[ 5] = (unsigned char)(width>> 8);
         infoHeader[ 6] = (unsigned char)(width>>16);
         infoHeader[ 7] = (unsigned char)(width>>24);
         infoHeader[ 8] = (unsigned char)(height    );
         infoHeader[ 9] = (unsigned char)(height>> 8);
         infoHeader[10] = (unsigned char)(height>>16);
         infoHeader[11] = (unsigned char)(height>>24);
         infoHeader[12] = (unsigned char)(1);
         infoHeader[14] = (unsigned char)(bytesPerPixel*8);
         
         return infoHeader;
      }


   private:
      const int height;
      const int width;
      std::vector<unsigned char> image;
   };
}; // namespace EVP


#endif // HEADER_IMAGE_BMP_HPP



