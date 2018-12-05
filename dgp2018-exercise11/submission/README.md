# Exercise 11 - Solving Laplace Equations II

Group 7: 
- Arne Sachtler
- Jannik Reichert
- Niklas Schmitz 

_(uniformly distributed workload)_

## Implementation Comments

1. The given framework computes the displacement for the vertices componentwise.
In our opinion this has the drawback, that the (very same) Laplacian of Laplacian matrix
has to be set up three times. We computed the update for all coordinates in the same
method and removed the auxiliary method.

2. We added a few buttons to the GUI to conveniently select the deformation parameters.
We can select if we want to use uniform or cotan weights for the deformations and whether
to use the thin-plate (L²x=b) or minimal surface (Lx=b) approach.

## Comparing Lx = b and L²x = b

The displacement can be interpreted as another mesh with the same graph structure but different vertex positions.
Therefore, we the use of the standard Laplacian imposes a minimal surface operation on the displacement mesh and the
squared Laplacian corresponds to the thin plate approach on the displacement mesh. Hence, for L the mean curvature of
the displacement mesh approximates zero and for L² the Gaussian curvature approximates zero. In consequence, we observe
a more global deformation for the L² approach, i. e. the displacement in a near locality of the fixed and shifted vertices
is smaller as for the L approach. However, the displacement variation is also smaller for the L² approach.

## Comparing Different Laplacian Weights

As in the previous assignments the uniform weighted Laplacian results in a bumpy mesh. We evaluated the difference on the
Max Planck mesh as the triangulation is highly irregular and we expect a larger difference between cotangent weighted and
uniform weighted Laplacian. We attached a screenshot of the result in "comparison_cotan_uniform.png". The deformation
based on the uniform Laplacian is worse as the approximation of the Laplacian is worse.

## Comparing with Physical Deformation

Collestion of our observations:
- The volume changes. This is not the case in the physical case, is it? It would lead to pressure differences.
- The surface never tears, no matter how much we deform.
- We can easily create self intersections, which do not occur in the real world.
- Deformation effects the complete mesh.

