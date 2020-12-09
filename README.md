# Distance Between Two Points On 3D Surface
In computer graphics, there are many known algorithms for defining the *distance between two points on a 3D surface*. So far, the most popular distance definitions are the **Geodesic distance** and the **Diffusion distance**, however, each of them has its advantages and drawbacks.

In general, the quality of a distance definition is measured by fulfillment of the following 8 properties:

1. **metric**: non-negative, satisfies the identity of indiscernibles, symmetric, and satisfies the triangle inequality;
2. **smooth**: smooth with respect to perturbations of $x$ and $y$, with no singularities except derivative discontinuity at $x$;
3. **locally isotropic**: approximately geodesic when y is near $x$;
4. **globally "shape-aware"**: reflects the overall shape of the surface when $y$ is far from $x$;
5. **isometry invariant**: does not change with isometric transformations of the surface; 
6. **insensitive to noise and topology**: does not change significantly with the addition of noise or changes to topology;
7. **practical to compute**: compute times between all pairs of points in common meshes take at most a few minutes; and
8. **parameter-free**: independent of any parameter that must be set differently for specific meshes or applications.

Unfortunately, none of the two popular methods satisfies all eight properties. But luckily, the **Biharmonic distance** was developed, which satisfies all eight properties.

|                                      | Geodesic distance | Diffusion distance | Biharmonic distance |
| ------------------------------------ | :---------------: | :----------------: | :-----------------: |
| 1) metric                            |         ✅         |         ✅          |          ✅          |
| 2) smooth                            |         ❌         |         ❌          |          ✅          |
| 3) locally isotropic                 |         ✅         |         ❓          |          ✅          |
| 4) globally "shape-aware"            |         ❌         |         ❓          |          ✅          |
| 5) isometry invariant                |         ✅         |         ✅          |          ✅          |
| 6) insensitive to noise and topology |         ❌         |         ✅          |          ✅          |
| 7) practical to compute              |         ✅         |         ✅          |          ✅          |
| 8) parameter-free                    |         ✅         |         ❌          |          ✅          |

✅ satisfied		❓not always satisfied		❌ not satisfied



### Biharmonic Distance

##### 1. Continuous definition

The **continuous biharmonic distance** is defined by 

$$d_{B}(x, y)^{2}=\sum_{k=1}^{\infty} \frac{\left(\phi_{k}(x)-\phi_{k}(y)\right)^{2}}{\lambda_{k}^{2}},$$

where $\phi_{k}(x)$ and $\lambda_{k}$ are the eigenfunctions and eigenvalues (resp.) of the positive definite Laplace-Beltrami operator

$$\Delta \phi_{k}(x)=\lambda_{k} \phi_{k}(x),$$

with $0=\lambda_{0}<\lambda_{1} \leq \lambda_{2} \ldots$



##### 2. Discrete definition

In the discrete case, the biharmonic distance is defined based on constructing a discrete Green’s function $g_d$ of the Bi-Laplacian using the formula

$$g_{B}(x, y)=\sum_{k=1}^{\infty} \frac{\phi_{k}(x) \phi_{k}(y)}{\lambda_{k}^{2}}.$$



Apply the common “cotangent formula” discretization of the Laplace-Beltrami differential operator on meshes, we get

$$\left(\Delta_{d} u\right)_{i}=\frac{1}{A_{i}} \sum_{j \in N e i(i)}\left(\cot \alpha_{i j}+\cot \beta_{i j}\right)\left(u_{i}-u_{j}\right),$$

where

* for a mesh function $u$ with $N$ vertices, $\left(\Delta_{d} u\right)_{i}$ denotes its discrete Laplacian evaluated at vertex $i$ (for $i = 1, 2, ..., N$)

* $A_i$ is the Voronoi area at the $i^{th}$ mesh vertex
* $α_{ij} , β_{ij}$ are the two angles supporting the edge connecting vertices $i$ and $j$

Let us denote the matrix of this linear transformation by $L_d$.



Let $A∈ \mathbb{R}^{N\times N}$ be the area/mass matrix of $u$, which is a diagonal matrix with $A_{ii} = A_{i}$.

Let $L_c$ be the conformal discrete Laplacian.

Then $L_{d}=A^{-1} L_{c}$.

The discrete Green’s function of the Bi-Laplacian (i.e. $g_d ∈ \mathbb{R}^{N\times N}$) can then be defined as the pseudo-inverse of $L_c A^{-1} L_c$.



Finally, the **discrete biharmonic distance** is defined by 

$$d_{B}\left(v_{i}, v_{j}\right)^{2}=g_{d}(i, i)+g_{d}(j, j)-2 g_{d}(i, j).$$



Libigl includes a wrapper for the discrete biharmonic distance$^{[1]}$ developed by Yilin Qu ([https://github.com/YILINQ/biharmonic_distance/blob/master/src/biharmonic_distance.cpp](https://github.com/YILINQ/biharmonic_distance/blob/master/src/biharmonic_distance.cpp)), exposing it through an Eigen-based API. The function

```C++
igl::biharmonic_distance(V,F,D);
```

computes the discrete biharmonic distances between each pair of vertices in V.

<img src="./image/cactus.png" alt="cactus" style="zoom:80%;" />

<img src="./image/lucy.png" alt="lucy" style="zoom:80%;" />

[Example 001](https://github.com/YILINQ/biharmonic_distance/blob/master/main.cpp) allows to interactively pick the source vertex and displays the distance using a periodic color pattern. (left) Biharmonic distance, (right) Geodesic distance.





### References

1. Y. Lipman, R. M. Rustamov, T. A. Funkhouser, [Biharmonic Distance](https://www.cs.princeton.edu/~funk/biharmonic.pdf), 2010.

