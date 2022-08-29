# MLOP - Manifold Reconstruction and Denoising from Scattered Data in High Dimension

This code implements the algorithm proposed in the paper
<ul>
 <li> 
[1] Faigenbaum-Golovin, S., & Levin, D. (2020). Manifold Reconstruction and Denoising from Scattered Data in High Dimension via a Generalization of $ L_1 $-Median. arXiv preprint arXiv:2012.12546.
  </li>
</ul>

If you find it useful, please cite our article.

In this paper, we present a method for denoising and reconstructing a low-dimensional
manifold in a high-dimensional space. We introduce a multidimensional extension of the
Locally Optimal Projection algorithm which was proposed by Lipman et al. in 2007 for
surface reconstruction in 3D. The high-dimensional generalization bypasses the curse of
dimensionality while reconstructing the manifold in high dimension. Given a noisy point-
cloud situated near a low dimensional manifold, the proposed solution distributes points
near the unknown manifold in a noise-free and quasi-uniformly manner, by leveraging a
generalization of the robust L1-median to higher dimensions. We prove that the non-
convex computational method converges to a local stationary solution with a bounded
linear rate of convergence if the starting point is close enough to the local minimum.
The effectiveness of our approach is demonstrated in various numerical experiments, by
considering different manifold topologies with various amounts of noise, including a case
of a manifold of different co-dimensions at different locations.


The files which should be in the directory:

manifold_reconstruction_MLOP.m    The main program. This program runs on the cylinder toy example. but other manifolds can be used also.
Data_Bunny.m                                                                
Data_Cone_Non_Uniform_Sample.m                                              
Data_Cylinder_Non_Uniform_Sample.m                                          
Data_Interface.m                                                            
Data_Orthogonal_Matrix_Uniform_Sample.m                                     
Data_RN_Cylinder_Non_Uniform_Sample.m                                       
Data_long_Cone_Non_Uniform_Sample.m                                         
Hd.m                                                                                                                                     
calculateFillDistancePerPoint.m                                             
calculateFillDistancePerSinglePoint.m                                       
calculateNorm.m                                                             
eta_function.m                                                                                                                   
main_manifold_reconstruction_MLOP_v1.m                                      
manifold_reconstruction_MLOP.m                                              
theta.m 

Troubleshooting:
- The problem should run as-is after everything is installed (creating "outdir" directory).
  If it doesn't run properly, there is probably a problem with the h1 or h2 not calibrated correctly. Please look on the med_data_min_points parameter, to see that there is optimal number of point from Q, and from P (i.e. not too many points, but also enough) - see the explanation for optimal h selection in the paper.
  
  
(C) Copyright: Shira Faigenbaum-Golovin

shirafaigen@gmail.com

https://www.tau.ac.il/~alecsan1/
