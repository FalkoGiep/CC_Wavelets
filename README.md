# CC_Wavelets
This is a python code to compute two-term as well as three-term connection coeffiicients for all wavelet families of finite support on a unbounded interval. This code is an fix and elaboration on an existing code by Manuels ( https://github.com/manuels/db_cc3 ), where three-term connection coefficients are calculated for Daubechies wavelets. Most of this code is based on the equations given in ( Besora, 2004 ), which in turn is based on the theory of ( Latto et al, 1992 ). 

## Theory
The theory is well outlined in ( Besora, 2004) or in the code of Manuels ( https://github.com/manuels/db_cc3 ). A summary is given here.
### Two-term connection coefficients
Say we have a scaling function ![equation](http://latex.codecogs.com/gif.latex?\phi) that corresponds to a filter h. the two-term connection coefficients are given by:

![equation](http://latex.codecogs.com/gif.latex?\Gamma_{i}^m=\int\phi(x)\frac{d^m\phi_i(x)}{dx^m}dx )

To calculate these connection coefficients we can use the scaling equation to end up with an eigen-value problem ( Goedecker, 2009 ):

![equation](http://latex.codecogs.com/gif.latex?\Gamma_{i}^m=2^m\sum_{j,k}h_jh_k\Gamma_{2i-j+k}^m )

Besides this a normalization equation is needed ( since an eigenvector is determined up to a constant). This is done by using the moment equation of the corresponding wavelet. In ( Goedecker, 2009 ) it is claimed that for all wavelets we would have the following:

![equation](http://latex.codecogs.com/gif.latex?\sum_ii^m\Gamma_i^m=m! )

This results in a over-determined system which we can solve with a least-squares algorithm.

### Three-term connection coefficients
the three-term connection coefficients are given by:

![equation](http://latex.codecogs.com/gif.latex?\Omega_{ij}^{mn}=\int\phi(x)\frac{d^m\phi_i(x)}{dx^m}\frac{d^n\phi_j(x)}{dx^n}dx )

Which results in a similair eigen-value problem by using the scaling equation ( Latto et al, 1992 ):

![equation](http://latex.codecogs.com/gif.latex?\Omega_{lm}^{mn}=2^{m+n-1}\sum_{i,j,k=0}^{N-1}h_ih_jh_k\Omega_{2l+j-i,2m+k-i}^{mn} )

It turns out this eigenvalue problem is rank-deficient ( Chen, 1996 ) and several extra homogenous equations need to be added. We can again use the moment equations:

![equation](http://latex.codecogs.com/gif.latex?0=\sum_jj^q\Omega_{j,1}^{mn} )

![equation](http://latex.codecogs.com/gif.latex?0=\sum_kk^p\Omega_{1,k}^{mn} )

Where ![equation](http://latex.codecogs.com/gif.latex?q<m,p<n ). Besides that we need a normalization equation:

![equation](http://latex.codecogs.com/gif.latex?m!n!=\sum_{j,k}j^mk^n\Omega_{j,k}^{mn} )

## Contents

### CC2.py
Script with the functions for calculating two-term connection coefficients.
### CC3.py
Script with the functions for calculating three-term connection coefficients.
### test.py
Script with the functions to validate the obtained three-term connection coefficients.
### cc3LIB.py
A database with some known three-term connection coefficients.


## Problems
If you have a problem you could do the following:
- Check if your filter is properly scaled. Some filters need to be multiplied/divided by ![equation](http://latex.codecogs.com/gif.latex?\sqrt{2} ).
- It is assumed that what is stated about the moments in ( Goedecker, 2009 ) is true, but no proof of this is given. It might be that is is not true for all wavelet families and you could another equation for the moments.





## References:
Besora, J., Galerkin Wavelet Method for Global Waves in 1D, 2004
Goedecker, S., Wavelets and their application for the solution of partial differential equations in physics, 2009
Latto, A., Resnikoff, H., Tenenbaum, E., The evaluation of connection coefficients of compactly supported wavelets, 1992
Chen, M., The Computation of Wavelet-Galerkin Approximation on a Bounded interval, 1996
