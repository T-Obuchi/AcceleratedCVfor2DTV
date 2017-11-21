# AcceleratedCVon2DTVLR
Approximate cross-validation for linear regression penalized by terms of *L1* and two-dimensional total variation.

This is free software, you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3 or above. See LICENSE.txt for details.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

# DESCRIPTION
Using the estimated explanatory variables *x* given the measument matrix *A* and the measurement result *y*, this program computes and returns an approximate leave-one-out error (LOOE) and its standard error for linear regression penalized by *L1* norm and two-dimensional total variation (TV).

# USAGE
```matlab
   [LOOE,ERR] = LOOEapprox_2DTV(x,y,A,Nx,Ny,lambda_T,delta,theta)
```
Inputs:
- *x*: Estimated explanatory variables (N=Nx*Ny dimensional vector). A two-dimensional (2D) image is expected in common cases.
- *y*: Measurement result (M dimensional vector)
- *A*: Measurement matrix (M*N dimensional matrix)
- *Nx*: One side length of *x* in 2D.
- *Ny*: Another side length of *x* in 2D.
- *lambda_T*: Regularization weight of TV
- *delta*: Softening constant of TV. Default value is 10^(-4).
- *theta*: Threshold to determine clusters induced by TV. Default value is 10^(-12).

Outputs:
- *LOOE*: Approximate value of the leave-one-out error
- *ERR*: Approximate standard error of the leave-one-out error

For more details, type help LOOEapprox_2DTV.

# REFERENCE
Tomoyuki Obuchi, Shiro Ikeda, Kazunori Akiyama, and Yoshiyuki Kabashima: "Accelerating cross-validation with total variation and its application to super-resolution imaging", arXiv: 1611.07197
