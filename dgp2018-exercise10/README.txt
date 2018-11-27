Exercise 10 - Parametrization and Minimal Surfaces

Group 7: Arne Sachtler, Jannik Reichert, Niklas Schmitz

Implementation Comments:

Exercise 1.2:
    We do not explicitely enforce the convexity constraints on the cotan weights, but rather multiply the texture coordinates
    by the unnormalized weights and divide the results by the sum of the weights. This way, we only need a single loop
    over the halfedges. Note that this implicitly enforces convexity of the weights used.
-- TODO ---
{ } warum 'almost closed' max so schlecht gemappt: boundary sehr klein (wenig vertices)
{ } stationaer bei implicit minimal surface nur bei uniform, da LGS gleich bleibt.
	nicht fuer cotan weights.
{ } was haben wir mit LX=0 right hand side gemacht (kurz erwaehnen)

