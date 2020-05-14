// Header only library to do union_find
// source:
// https://www.geeksforgeeks.org/union-find-algorithm-set-2-union-by-rank/

#ifndef HEADER_UNION_FIND_HPP
#define HEADER_UNION_FIND_HPP
#include <iostream>
#include <vector>

namespace evp {

class Union_find {
public:
   Union_find(const int n) : n(n) {
      subsets.resize(n);
      for (int v = 0; v < n; ++v) {
         subsets[v].parent = v;
         subsets[v].rank = 0;
      }
   }

   // A utility function to find set of an element i
   // (uses path compression technique)
   int find(int i) {
      // find root and make root as parent of i (path compression)
      if(subsets[i].parent != i) {
         subsets[i].parent = find(subsets[i].parent);
      }
   
      return subsets[i].parent;
   }
   
   // A function that does union of two sets of x and y
   // (uses union by rank)
   void Union(int x, int y) {
      int xroot = find(x);
      int yroot = find(y);
   
      // Attach smaller rank tree under root of high rank tree
      // (Union by Rank)
      if (subsets[xroot].rank < subsets[yroot].rank) {
         subsets[xroot].parent = yroot;
      } else if (subsets[xroot].rank > subsets[yroot].rank) {
          subsets[yroot].parent = xroot;
      } else {
         // If ranks are same, then make one as root and increment
         // its rank by one
         subsets[yroot].parent = xroot;
         subsets[xroot].rank++;
      }
   } 
private:

   struct subset { 
      int parent; 
      int rank; 
   };

   std::vector<subset> subsets;

   const int n;
};


   
};// end namespace evp


#endif // HEADER_UNION_FIND_HPP
