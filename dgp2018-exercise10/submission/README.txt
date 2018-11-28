Exercise 10 - Parametrization and Minimal Surfaces

Group 7: Arne Sachtler, Jannik Reichert, Niklas Schmitz (uniformly distributed workload)

Implementation Comments:

Exercise 1.2:
    We do not explicitly enforce the convexity constraints on the cotan weights, but rather multiply the texture coordinates
    by the unnormalized weights and divide the results by the sum of the weights. This way, we only need a single loop
    over the halfedges. Note that this implicitly enforces convexity of the weights used.


Results:

Poor Texture on Max Head:

We observe that the mapping of the head model the unit circle results in a bad textured result. Looking at the 2D parameter
domain reveals that the entire head is maps very close to the center of the circle. This is because we have a small boundary
consisting of a small amount of vertices which are connected to each other.

Comparison of Uniform and Cotangent Weighted Minimal Surfaces:


-- TODO ---
{x} warum 'almost closed' max so schlecht gemappt: boundary sehr klein (wenig vertices)
{ } stationaer bei implicit minimal surface nur bei uniform, da LGS gleich bleibt.
	nicht fuer cotan weights.
{ } was haben wir mit LX=0 right hand side gemacht (kurz erwaehnen)

