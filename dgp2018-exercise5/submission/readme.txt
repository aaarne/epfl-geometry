Exercise 5: Surface Normals

The computation of normals, curvatures and the laplacian are implemented by traversing the mesh data structure.
Instead of using the circulators directly in the c-style version, we decided to use them in more readably for-each
expressions. Note that the for each loop internally does the same thing as creating a do-while loop and incrementing
the pointer manually would do.

For the small sphere we get:
    Min Uniform Laplace value is: 0.0101589
    Max Uniform Laplace value is: 0.0168597
    Min Laplace-Beltrami curvature value is: 0.949649
    Max Laplace-Beltrami curvature value is: 1.08865
    Min Gauss curvature value is: 0.91041
    Max Gauss curvature value is: 1.04125

For the large sphere we get:
    Min Uniform Laplace value is: 0.10159
    Max Uniform Laplace value is: 0.168598
    Min Laplace-Beltrami curvature value is: 0.0949652
    Max Laplace-Beltrami curvature value is: 0.108865
    Min Gauss curvature value is: 0.0091041
    Max Gauss curvature value is: 0.0104125

Comparison of the Uniform Laplacian:
As expected, the uniform Laplacian increases by a factor of 10 as the sphere radius increases by a factor of ten.
Therefore, the distance vectors of each vertex to the surrounding vertices are scaled by the factor of 10. As the uniform
Laplacian is an average of the outgoing halfedges of a vertex and it our resizing of the sphere we scale every vector by
a factor of 10, these halfedges get also scaled by 10. In the average we can pull the constant factor out of the sum.
Consequently, the uniform Laplacian is scaled by the same factor.

Comparison of the Laplace-Beltrami Curvature:
First, note that the angles in a mesh are invariant w.r.t. scaling. Again, a factor of 10 appears in the sum, but we also
have a factor of 10^2=100 in the denominator. Therefore, the total value for the Laplace-Beltrami curvature is scaled
by a factor of one tenth.

Contributions:
- Jannik: 1/pi
- Niklas: 34.08%
- Arne: 34.08%
