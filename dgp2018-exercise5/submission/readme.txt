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
