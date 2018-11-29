Exercise 10 - Parametrization and Minimal Surfaces

Group 7: Arne Sachtler, Jannik Reichert, Niklas Schmitz (uniformly distributed workload)

Implementation Comments:

Exercise 1.2:
    We do not explicitly enforce the convex combination constraints on the cotan weights, 
    but rather multiply the texture coordinates by the unnormalized weights and 
    divide the results by the sum of the weights. Thus, we only need a single loop
    over the halfedges. Note that this implicitly enforces a convex combination
    of the weights used. We assert non-negativity for all weights explicitly and observe
    that it does not fail, which is also obvious looking at the code of calc_edge_weights().

Exercise 1.3:
    We tried both approaches:
        1. Solving the simple (n+m)x(n+m) having the one-hot encoded rows for boundary vertices and
        2. Solving the reduced nxn system as provided in the slides.
    We verified successfully, that the reduced system (2.) and the complete system (1.) yield the same results.
    However, comparing the run-time shows that the overall runtime of the reduced system is slower as we need
    more time to set up the matrices and the need to perform a mapping from vertex to matrix indicies.
    On the Max Planck mesh the full system (1.) need ~85ms in release profile, while the reduced
    system (2.) needs 110ms. So even though the system to be solved is smaller in terms of dimensions the
    runtime is longer. That's why we keep the (n+m)x(n+m) approach for the final solution as it is (surprisingly)
    faster and the code is simpler and therefore more beautiful.

Exercise 2:
	We adapt the homogeneous linear equation system of the form LX=0 by adding the boundary
	constrains to the right-hand side of the system. Additionally, we set the rows of the L matrix
	for the boundary vertices to 1 at a single point.

Results:

Poor Texture on Max Head:

	We observe that the mapping of the head model onto the unit circle results in poor texture quality.
	Looking at the 2D parameter domain reveals that the most of the head maps very close to the 
	center of the circle. This is because we have a comparatively small boundary
	consisting of a small amount of vertices which are connected to each other.

Comparison of Uniform and Cotangent Weighted Minimal Surfaces:

	The discrete Laplacian matrix is only depending on vertex connectivity if we use uniform weights.
	Consequently, the solution of the linear equation system is always the same if we solve for the 
	minimal surface on the mesh multiple times.

	In the folder "minimal_uniform" you can see the initial state of cylinder3 and the result of the minimal
	surface operation on the mesh after the first and the second run. We observe that indeed the solution of
	the minimal of the initial surface is always the same no matter how often we execute the procedure.

	On the other hand, the folder "minimal_cotan" shows the results of applying the minimal surface opration multiple
	times using cotan weights. We observe that the solution is not stationary as the system of equations depends not
	only on the graph structure of the mesh, but also on the geometry. Additionally, after the 4th application a strange
	singularity appears, where the cylinder has zero diameter everywhere except at the boundary circles.
	Iterating even longer results in a periodic behaviour yielding the results shown in "after_n.png" and "after_n+1.png"
	alternatingly.
	Also looking at the amount of edges having 0 weights shows an  characteristic behaviour using cotan weights when
	applying the minimal surface operator multiple times on the cylinder3 model:
        0 out of 825 cotan weights were 0. (1. iteration)
        0 out of 825 cotan weights were 0.
        30 out of 825 cotan weights were 0.
        30 out of 825 cotan weights were 0.
        58 out of 825 cotan weights were 0.
        60 out of 825 cotan weights were 0.
        53 out of 825 cotan weights were 0.
        45 out of 825 cotan weights were 0.
        30 out of 825 cotan weights were 0.
        105 out of 825 cotan weights were 0.
        90 out of 825 cotan weights were 0. (11. iteration)
        90 out of 825 cotan weights were 0.
        90 out of 825 cotan weights were 0. ... (goes on with 90 forever)

