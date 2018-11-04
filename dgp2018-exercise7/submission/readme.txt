Exercise 7 - Smoothing
----------------------

General: We documented the approaches how the implementation assignments are solved in the code.
For exercise 2 (implicit smoothing) we use the following equation to compute the updated vertex positions:

 /  -1                          \  (t+1)   -1      (t)
 | D   - delta * t * lambda * M | P     = D    *  P
 \                              /

 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  ^^^^^  ^^^^^^^^^^^^
                A                  X    =     B

where we create the sparse matrix A and the dense matrix B.


Questions
---------

Do you experience analogous behavior for surfaces regarding the Gage-Hamilton-Grayson theorem?

TODO


Feature enhance on bad and nice face mesh

We observe, that the feature enhancing of the bad mesh fails and completely distorts the resulting output mesh.
The cotan feature enhancement performs better, but only produces a slight feature enhancement.
For the nice mesh both algorithms result in similar results.
Screenshots of the results:
- bad_mesh_uniform_enhancement.png: Bad mesh using uniform feature enhancement
- bad_mesh_cotan_enhancement.png: Bad mesh using cotan feature enhancement
- nice_mesh_uniform_enhancement.png: Nicd mesh using uniform feature enhancement
- nice_mesh_cotan_enhancement.png: Nice mesh using cotan feature enhancement

Feature enhancement results:

- large_ear_bunny.png: We created and interesting new creature that looks a bit like a sharp-corner bunny and has large ears
- spooky_max.png: The halloween version of Max Planck
- large_nose_fox.png: A feature-enhanced version of our own scanned model
- UNIL_fox_original.png: The original version of the UNIL fox (for reference)
