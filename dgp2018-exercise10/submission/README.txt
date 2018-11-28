Exercise 10 - Parametrization and Minimal Surfaces

Group 7: Arne Sachtler, Jannik Reichert, Niklas Schmitz (uniformly distributed workload)

Implementation Comments:

Exercise 1.2:
    We do not explicitly enforce the convex combination constraints on the cotan weights, 
    but rather multiply the texture coordinates by the unnormalized weights and 
    divide the results by the sum of the weights. Thus, we only need a single loop
    over the halfedges. Note that this implicitly enforces a convex combination
    of the weights used. We assert non-negativity for all weights explicitly.

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
	TODO: add fancy plot of percentage of zero weight for cotan



-- TODO ---
{x} warum 'almost closed' max so schlecht gemappt: boundary sehr klein (wenig vertices)
{ } stationaer bei implicit minimal surface nur bei uniform, da LGS gleich bleibt.
	nicht fuer cotan weights.
{ } was haben wir mit LX=0 right hand side gemacht (kurz erwaehnen)
