Exercise 8: Remeshing

Group 7: Arne Sachtler, Jannik Reichert, Niklas Schmitz

As usual, codes are implemented as described on the exercise sheet, details are in the comments.
Note:
- we added a new function in mesh_processing.cpp: calc_max_curvature()
	- called in compute_mesh_properties, it is stored as a vertex property
	



===================================================================
- kommentare erklaeren code
- note that moved computation of the curvature to own method in order to debug the max curvature
- When calculating the laplacian projection onto the tangent plane, we use an Eigen hyperplane to improve readability
- added shader for target length and max curvature
- if max curvature calculation fails with neg sqrt, we set it to H

Average Remeshing:
Q:
Does this algorithm lead to a stationary solution after finite number of steps, i.e. does the mesh not change anymore when more remeshing steps are applied?

Short Answer: No.

Long Answer: If we try many iteration of the remeshing loop (>1000) the mesh does not converge to a stationary distribution. In particular, the algorithm keeps splitting,
collapsing and flipping edges. We count the number of edges treated in each iteration and observe that is does not tends towards zero. More accurately, we observe that
even after the 1000th iteration still ~3000 edges are split and collapsed and ~60 edges are flipped in total per iteration. We have several arguments explaining the 
non-convergence:
1. Consider the splitting and collapsing of edges. Long edges are split if their length is above four thirds the target length and short edges are collapsed if below four
fifth the target length. Let us show the non-convergence with an example: Consider the target length being 1 and two vertices connected by an edge of length 1.5.
The edge will be split as 1.5 > 1.333 resulting in a new vertex and two edges of length 0.75. As 0.75 < 0.8 = 4/5 the new vertex will be collapsed and we get back
the edge of length 1.5. This process continues forever and the solution is therefore periodic and non-convergent. Note that this argument assumes the non-existence of 
edge flip and tangential relaxation, but it gives a first intuition of why the solution is non-stationary.
2. We do not change the boundary vertices. Note that all operations except the collapsing keep the sum of valence constant. Now also note that the collapsing action 
reduces the sum of valence by a multiple of two. We therefore cannot converge to a point of equally distributed valence (as in the football case) and the edge flip 
operation will move point of high or low valence over the mesh without convergence.





