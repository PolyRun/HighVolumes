// Header only for split

#include <string>
#include <vector>
#include <sstream>

#ifndef HEADER_SPLIT_HPP
#define HEADER_SPLIT_HPP

std::vector<std::string> split(std::string &s, char separator =',') {
    std::vector<std::string> res;
    std::istringstream f(s);
    std::string tmp;
    while (getline(f, tmp, separator)) {
        res.push_back(tmp);
    }
    return res;
}

#endif // HEADER_SPLIT_HPP
