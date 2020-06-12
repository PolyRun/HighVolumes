# Papers

## Algorithm

The specific algorithm we would like implement to is "Algorithm 1" given in this paper:

[A Fast and Practical Method to Estimate Volumes of Convex Polytopes](https://arxiv.org/abs/1401.0120)
(Cunjing Ge, Feifei Ma, Jian Zhang 2013)


## Survey
[How to compute the volume in high dimension?](https://www.math.tamu.edu/~rojas/simonvitzvolumehigh.pdf)
(Miklos Simonovits 2003)

## Find Ellipsoids

In Chapter 3 of this book, they discuss the Ellipsoid method to check the feasibility for linear inequalities.

Theorem 3.1.9 states that there is exactly one Ellipsoid E of minimal volume for any convex body K, such that E contains K (minimal outer ellipsoid).

Given various separation oracles for K, we can either use a `deep cut` to find the minimal ellipsoid faster, or a `shallow cut` for weaker separation oracles and find the ellipsoid a bit slower, though still in polynomial time.

[Geometric Algorithms and Combinatorial Optimization](https://www.springer.com/gp/book/9783642782428)
(Grötschel, Martin, Lovasz, Laszlo, Schrijver, Alexander 1993)
