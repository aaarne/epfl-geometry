# Exercise 11 - Solving Laplace Equations II

Group 7: 
- Arne Sachtler
- Jannik Reichert
- Niklas Schmitz 

_(uniformly distributed workload)_

## Implementation Comments

1. The given framework computes the displacement for the vertices component-wisely.
In our opinion this has the drawback, that the (very same) Laplacian of Laplacian matrix
has to be set up three times. We computed the update for all coordinates in the same
method and removed the auxiliary method.

2. We added a few buttons to the GUI to conveniently select the deformation parameters.
We can select if we want to use uniform or cotan weights for the deformations and whether
to use the thin-plate (L²x=b) or minimal surface (Lx=b) approach.

## Comparing Lx = b and L²x = b

TODO

## Comparing Different Laplacian Weights

TODO

## Comparing with Physical Deformation

TODO, some ideas:
- The volume changes. This is not the case in the physical case, is it? It would lead to pressure differences.
- The surface never tears, no matter how much we deform.
- We can easily create self intersections, which do not occur in the real world.
- Deformation effects the complete mesh. If we fix the tail and deform an ear of the bunny, the complete mesh (except the tail)
changes. I believe in a physical world it is possible to only deform the ear of a bunny.

