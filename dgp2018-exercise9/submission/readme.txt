Exercise 9: Solving Laplace Equations 1

Group 7: Arne Sachtler, Jannik Reichert, Niklas Schmitz
uniformly distributed workload

Solving Laplace Equations
-------------------------

We construct the discrete Laplacian matrix using the formulas given in the lecture. Recall the factorization of the
laplacian L=DM. As we solve the system of equation to zero (we want a harmonic function which requires the laplacian
to be zero at every point), we see that the diagonal matrix D does not influence the solution. Therefore, we compute
only the matrix M using the following equations if there is no constaint on the vertex i:

m_ij = cotan[e_ij] if j ∈ N1(i),
m_ii = Σcotan[i_ik] for k ∈ N1(i),
and 0 otherwise.

If there is a constraint on the vertex i we set

m_ii = 1
and m_ij = 0 ∀j

futher set the i_th element of rhs to 1 if i is the second vertex with constraint. All other elements of rhs remain zero.

Generating Edges at Isolines of Harmonic Functions
--------------------------------------------------

Given the isovalues on the interval borders a and b, we first compute the minimum and maximum value of the isolines.
Then, the lowest index of the isoline passing through the interval ab is

             /  min(a,b) - l \                           /  max(a,b) - l \
first = ceil | ------------- |          and last = floor | ------------- |
             \ interval_size /                           \ interval_size /

For the coordinates of the isolines we use interpolation on the vector spanned by two vertices of a triangle, respectively.
We compute the interpolation ratio by the isovalue of the isoline and the isovalues on the two boundaries of the interval.

Results
-------

We tried evaluated the implementation with three different experiments, which you can see in the attached screenshots.

- octahedron_corner_constraints.png: We set the four constrained points on four corners of the octahedron such that 1/0
and 2/3 are on opposite corners. We can observe nice level sets.

- max_boundary_contraintes.png: We set the four constraints on the boundary of the Max Planck mesh.

- isolines_on_rabbit.png: We set a constraint of the tail and nose of the bunny and two constraints on the ear tips.
