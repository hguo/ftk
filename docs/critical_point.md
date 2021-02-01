# Tutorial: Critical Point Tracking with FTK

Critical points---vanishing locations in a vector field---are the most important family of features in scientific data.  Critical points are the key constitiuents of vector field topology and essentially determine the characteristics of flow transport such as sources, sinks, and saddles.   Critical points are also in the gradient field of scalar functions, representing minimum, maximum, and saddles in a scalar field.  Read [this preprint](https://arxiv.org/abs/2011.08697) for more details on critical point tracking in FTK. 

## Critical point extraction

We use a numerical method to locate critical points in the input vector field data.  Without the loss of generality, we describe the case of 2D vector fields and assume the vector field is *piecewise linear* (PL), which implies that the data defined on a simplicial (triangular) mesh.  We will discuss how our method adapt to non-simplicial meshes in following sections.  

To extract critical points in a 2D PL vector field, one can iterate each triangular cell and test if the cell encircles a critical point.  

To be continued...
