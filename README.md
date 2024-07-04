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
