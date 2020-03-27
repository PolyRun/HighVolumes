

/**
 * \brief transform input polytope s.t. it fits in unit ball and has 1/2n * unit ball inside it
 * \param P: input polytope
 * \param Q: output polytope
 * \param det: output: the determinant of the linear transformation applied to P in order to get Q
 **/
void preprocess(Polytope *P, Polytope **Q, double *det);
