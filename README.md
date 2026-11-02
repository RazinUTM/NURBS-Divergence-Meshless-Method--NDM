# NURBS-Divergence-Meshless-Method--NDM

N. Rasin a,e  
H. Hirol a,e  
A. R. Zainal Abidin b, e, *  
M. H. Mokhtaram c,e  
M. A. Mohd Nor c,e  
A. Y. Mohd Yassin d,e  

a) Malaysia-Japan International Institute of Technology (MJIIT), Universiti Teknologi Malaysia, Kuala Lumpur, Malaysia  
b) Faculty of Civil Engineering, Universiti Teknologi Malaysia, Johor Bahru, Malaysia  
c) Faculty of Engineering and Life Sciences, Universiti Selangor, Malaysia  
d) School of Energy, Geoscience, Infrastructure and Society (EGIS), Heriot-Watt University Malaysia  
e) Malaysian Society for Numerical Methods (MSNM)  
* Corresponding author (A.R. Zainal Abidin, arazin@utm.my)  

Executive Summary:  

Inspired by the earlier works on MFV and drawing from our experience in MM and coupling NURBS with 
radial point interpolation method (RPIM), this work proposes an extension of the divergence theorem 
to meshless methods, with geometry represented by NURBS. As will be shown, the use of NURBS not only 
ensures geometric exactness but also facilitates the treatment of the product of the normal direction and
spatial differential terms, due to the rectangular nature of the parameter space. The proposed method, 
termed the NURBS-Divergence-Meshless (NDM) method, differs from previous MFVs in the following aspects:

i. Utilization of NURBS for both geometric exactness and simplified handling of normal directions and spatial differential terms.  
ii. Implementation of Gauss numerical integration beyond the one-point Gaussian rule.  
iii. Insensitivity to overlapping integration cells.  

Here we present 3 examples of it application:  
Ex1) Fan Blade Problem  
Ex2) Quarter-annular Disc  
Ex3) Cellular Beam  

## Code authorship

The MATLAB code in this repository was written by Ahmad Razin Zainal Abidin (Universiti Teknologi Malaysia) as part of the Malaysian Society for Numerical Methods (MSNM) work on NDM. The method itself is the joint work of the six authors listed above.

## How to cite

If you use this code, please cite the article:

N. Rasin, H. Hirol, A. R. Zainal Abidin, M. H. Mokhtaram, M. A. Mohd Noor and A. Y. Mohd Yassin (2025). NURBS-Divergence-Meshless (NDM) method. *Computers & Mathematics with Applications*, 196, 151–171. https://doi.org/10.1016/j.camwa.2025.07.011

A machine-readable citation is in [`CITATION.cff`](CITATION.cff).

## Licence

MIT Licence, see [`LICENSE`](LICENSE). Copyright (c) 2021–2024 Ahmad Razin Zainal Abidin and Universiti Teknologi Malaysia.

## Disclaimer

Research code, provided as is and without warranty. Results must be checked independently before they are used in design or safety decisions.

## AI-assisted maintenance

The 2026 update that added the licence, the citation file and the copyright headers was prepared with AI assistance and reviewed by the author. No numerical code was changed. See [`AI_USAGE.md`](AI_USAGE.md).
