Exercise 8: Remeshing

Group 7: Arne Sachtler, Jannik Reichert, Niklas Schmitz

As usual, codes are implemented as described on the exercise sheet, details are in the comments.
A few notes:

Tangential Smoothing
- 

Adaptive Remeshing
- we added a new function in mesh_processing.cpp: calc_max_curvature()
	- called in compute_mesh_properties, it is stored as a vertex property
	- we encountered rare (5 out of 10000 vertices on Max Planck) problems with negative 
	  squareroots when calculating the max curvature from mean curvature H and gaussian 
	  curvature K, i.e. in the formula
	  k_max = H + sqrt(H^2 - K)
	  In this case we proceed by setting it just to the mean_curvature, i.e. k_max = H


===================================================================
- kommentare erklaeren code
- note that moved computation of the curvature to own method in order to debug the max curvature
- When calculating the laplacian projection onto the tangent plane, we use an Eigen hyperplane to improve readability
- added shader for target length and max curvature
- if max curvature calculation fails with neg sqrt, we set it to H

Average Remeshing:
Q:
Does this algorithm lead to a stationary solution after finite number of steps, i.e. does the mesh not change anymore when more remeshing steps are applied?
A:
No





