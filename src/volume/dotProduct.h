// imported by volume.h

// Reference impl:
FT dotProduct_ref(const FT* u, const FT* v, const int n);

// First vectorization:
FT dotProduct_vec1(const FT* u, const FT* v, const int n);
// comment:
//  seems to be slower than ref, takes more than 2x time

FT dotProduct_auto1(const FT* u, const FT* v, const int n);
