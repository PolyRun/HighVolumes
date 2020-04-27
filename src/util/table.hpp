// Header only library to print a file to std::cout

#ifndef HEADER_TABLE_HPP
#define HEADER_TABLE_HPP
#include <iostream>
#include <map>
#include <vector>

class Table {
public:
   Table(const std::vector<std::string> &keys)
   : keys(keys) {
      col_size.resize(keys.size());
      for(int i=0;i<keys.size();i++) {
         column[keys[i]] = i;
	 col_size[i] = keys[i].size();
      }
      rows.push_back(keys);
   }

   void addRow(const std::map<std::string, std::string> &row) {
      std::vector<std::string> next;
      for(int i=0;i<keys.size();i++) {
         auto it = row.find(keys[i]);
	 if(it == row.end()) {
	    next.push_back("-");
	    col_size[i] = std::max(col_size[i],1);
	 }else{
	    next.push_back(it->second);
	    col_size[i] = std::max(col_size[i],(int)(it->second.size()));
	 }
      }
      rows.push_back(next);
   }

   void print() {
      for(int r=0;r<rows.size();r++) {
         for(int i=0;i<keys.size();i++) {
	    const std::string &v = rows[r][i];
	    std::cout << v;
	    for(int j=0;j<col_size[i]-v.size()+1;j++) {std::cout << " ";}
	 }
	 std::cout << "\n";
      }
   }
private:
   std::vector<std::string> keys;
   std::map<std::string,int> column;
   std::vector<int> col_size;
   std::vector<std::vector<std::string>> rows;
};


#endif // HEADER_TABLE_HPP
