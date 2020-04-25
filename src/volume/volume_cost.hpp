
#ifndef HEADER_VOLUME_COST_HPP
#define HEADER_VOLUME_COST_HPP

#include "../util/performance_counter.hpp"

typedef void (*xyz_cost_f) (const int);
void xyz_f1_cost(const int n);

// ----------------------------------------------- linalg
typedef void (*dotProduct_cost_f)(const int);
void dotProduct_cost_ref(const int n);
typedef void (*vectorNorm_cost_f)(const int);
void vectorNorm_cost_ref(const int n);

typedef void (*Ball_intersectCoord_cost_f)(const int);
void Ball_intersectCoord_cost_ref(const int n);

// ----------------------------------------------- volume

typedef void (*volume_cost_f)(const int, const int, const void**,const Body_T**);
void volume_cost_ref(const int n, const int bcount, const void** body, const Body_T** type);

typedef void (*walk_cost_f)(const int, int bcount, const void**, const Body_T**);
void walk_cost_ref(const int n, int bcount, const void** body, const Body_T** type);
void walkCoord_cost_ref(const int n, int bcount, const void** body, const Body_T** type);

#endif // HEADER_VOLUME_COST_HPP
