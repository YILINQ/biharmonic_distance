# Distance Between Two Points On 3D Surface
In computer graphics, there are many known algorithms for defining the *distance between two points on a 3D surface*. So far, the most popular distance definitions are the **Geodesic distance** and the **Diffusion distance**, however, each of them has its advantages and drawbacks.

In general, the quality of a distance definition is measured by fulfillment of the following 8 properties:

1. **metric**: non-negative, satisfies the identity of indiscernibles, symmetric, and satisfies the triangle inequality;
2. **smooth**: smooth with respect to perturbations of <img src="./svgs/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode" align=middle width=9.3949878pt height=14.1552444pt/> and <img src="./svgs/deceeaf6940a8c7a5a02373728002b0f.svg?invert_in_darkmode" align=middle width=8.64922575pt height=14.1552444pt/>, with no singularities except derivative discontinuity at <img src="./svgs/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode" align=middle width=9.3949878pt height=14.1552444pt/>;
3. **locally isotropic**: approximately geodesic when y is near <img src="./svgs/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode" align=middle width=9.3949878pt height=14.1552444pt/>;
4. **globally "shape-aware"**: reflects the overall shape of the surface when <img src="./svgs/deceeaf6940a8c7a5a02373728002b0f.svg?invert_in_darkmode" align=middle width=8.64922575pt height=14.1552444pt/> is far from <img src="./svgs/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode" align=middle width=9.3949878pt height=14.1552444pt/>;
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

<p align="center"><img src="./svgs/05d549fdd853244cfb6b0ad654f0d04d.svg?invert_in_darkmode" align=middle width=242.0901483pt height=46.5728307pt/></p>

where <img src="./svgs/b4044b42d577e0cbae46cc7ca353197d.svg?invert_in_darkmode" align=middle width=40.06290915pt height=24.657534pt/> and <img src="./svgs/cf39565086d308d92ed10730aba2a5bf.svg?invert_in_darkmode" align=middle width=16.855113pt height=22.8310566pt/> are the eigenfunctions and eigenvalues (resp.) of the positive definite Laplace-Beltrami operator

<p align="center"><img src="./svgs/1fb649907a282d20464355469371fa66.svg?invert_in_darkmode" align=middle width=137.9853717pt height=16.438356pt/></p>

with <img src="./svgs/f41dbd43745211255cdd5848f4d4edaf.svg?invert_in_darkmode" align=middle width=146.78037825pt height=22.8310566pt/>



##### 2. Discrete definition

In the discrete case, the biharmonic distance is defined based on constructing a discrete Green’s function <img src="./svgs/3f5518c046134f19d4b3f895f8a89a18.svg?invert_in_darkmode" align=middle width=14.6836602pt height=14.1552444pt/> of the Bi-Laplacian using the formula

<p align="center"><img src="./svgs/c64f08be4a156f02df104725b7754554.svg?invert_in_darkmode" align=middle width=193.74923535pt height=45.2741091pt/></p>



Apply the common “cotangent formula” discretization of the Laplace-Beltrami differential operator on meshes, we get

<p align="center"><img src="./svgs/b381ccfb9aa3ba55ed611bb32f307aab.svg?invert_in_darkmode" align=middle width=349.89933165pt height=45.00203565pt/></p>

where

* for a mesh function <img src="./svgs/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.4102734pt height=14.1552444pt/> with <img src="./svgs/f9c4988898e7f532b9f826a75014ed3c.svg?invert_in_darkmode" align=middle width=14.99998995pt height=22.4657235pt/> vertices, <img src="./svgs/5342555c02fffd4f41d5358a95d1adf3.svg?invert_in_darkmode" align=middle width=48.210261pt height=24.657534pt/> denotes its discrete Laplacian evaluated at vertex <img src="./svgs/77a3b857d53fb44e33b53e4c8b68351a.svg?invert_in_darkmode" align=middle width=5.6632257pt height=21.6830097pt/> (for <img src="./svgs/e55cab2d4adebe22b11cffc3e2ba805d.svg?invert_in_darkmode" align=middle width=94.6355883pt height=22.4657235pt/>)

* <img src="./svgs/4ebf880807deff5796460f39aea46f80.svg?invert_in_darkmode" align=middle width=16.9796979pt height=22.4657235pt/> is the Voronoi area at the <img src="./svgs/3def24cf259215eefdd43e76525fb473.svg?invert_in_darkmode" align=middle width=18.3250452pt height=27.9124395pt/> mesh vertex
* <img src="./svgs/17aabd8b859668124255fdc0fb94421b.svg?invert_in_darkmode" align=middle width=49.45218465pt height=22.8310566pt/> are the two angles supporting the edge connecting vertices <img src="./svgs/77a3b857d53fb44e33b53e4c8b68351a.svg?invert_in_darkmode" align=middle width=5.6632257pt height=21.6830097pt/> and <img src="./svgs/36b5afebdba34564d884d347484ac0c7.svg?invert_in_darkmode" align=middle width=7.710417pt height=21.6830097pt/>

Let us denote the matrix of this linear transformation by <img src="./svgs/c930bbff186713a73e826098c5d1b889.svg?invert_in_darkmode" align=middle width=18.03032055pt height=22.4657235pt/>.



Let <img src="./svgs/88cceb4989ce7185889653ac6717d534.svg?invert_in_darkmode" align=middle width=77.8584378pt height=27.6567522pt/> be the area/mass matrix of <img src="./svgs/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.4102734pt height=14.1552444pt/>, which is a diagonal matrix with <img src="./svgs/2be14126b642e853922b6cff93d8e3c4.svg?invert_in_darkmode" align=middle width=61.3498017pt height=22.4657235pt/>.

Let <img src="./svgs/5c14406eca8c4c05451828ed352eb491.svg?invert_in_darkmode" align=middle width=17.0618943pt height=22.4657235pt/> be the conformal discrete Laplacian.

Then <img src="./svgs/80e44676e0de1c28f8ddc1d9d25528a3.svg?invert_in_darkmode" align=middle width=87.80901855pt height=26.7617526pt/>.

The discrete Green’s function of the Bi-Laplacian (i.e. <img src="./svgs/008f0ab09a0630e0d9b188daf71a67fe.svg?invert_in_darkmode" align=middle width=81.03521745pt height=27.6567522pt/>) can then be defined as the pseudo-inverse of <img src="./svgs/3ccc0649de7307c25d7cc87a073fef7b.svg?invert_in_darkmode" align=middle width=64.92296085pt height=26.7617526pt/>.



Finally, the **discrete biharmonic distance** is defined by 

<p align="center"><img src="./svgs/ec779b818a88aa9e1098200222f50e51.svg?invert_in_darkmode" align=middle width=299.29501965pt height=20.3835786pt/></p>



Libigl includes a wrapper for the discrete biharmonic distance<img src="./svgs/ecb2edc691a4d089460ad80188b5da8f.svg?invert_in_darkmode" align=middle width=13.9955178pt height=29.190975pt/> developed by Yilin Qu ([https://github.com/YILINQ/biharmonic_distance/blob/master/src/biharmonic_distance.cpp](https://github.com/YILINQ/biharmonic_distance/blob/master/src/biharmonic_distance.cpp)), exposing it through an Eigen-based API. The function

```C++
igl::biharmonic_distance(V,F,D);
```

computes the discrete biharmonic distances between each pair of vertices in V.

<img src="./image/cactus.png" alt="cactus" style="zoom:80%;" />

<img src="./image/lucy.png" alt="lucy" style="zoom:80%;" />

[Example 001](https://github.com/YILINQ/biharmonic_distance/blob/master/main.cpp) allows to interactively pick the source vertex and displays the distance using a periodic color pattern. (left) Biharmonic distance, (right) Geodesic distance.



##### 3. Approximate discrete computation

Under the consideration of runtime, we may approximate <img src="./svgs/d68565aaa89193af9c2a54713859679b.svg?invert_in_darkmode" align=middle width=58.0059876pt height=24.657534pt/> by computing

<p align="center"><img src="./svgs/101979c7c8df113a315890e3d6369c59.svg?invert_in_darkmode" align=middle width=242.09012685pt height=48.18280005pt/></p>

That is to compute the first <img src="./svgs/d6328eaebbcd5c358f426dbea4bdbf70.svg?invert_in_darkmode" align=middle width=15.13700595pt height=22.4657235pt/> eigenvectors of the discrete Laplacian <img src="./svgs/8c5dbfcdf090de6e4792d3fde07eb146.svg?invert_in_darkmode" align=middle width=95.9013858pt height=22.8310566pt/>, which boils down to solving the generalized eigenvalue problem

<p align="center"><img src="./svgs/e51cf2c5dfe176066306391a120c0a36.svg?invert_in_darkmode" align=middle width=110.1384669pt height=14.6118786pt/></p>



Libigl includes a wrapper for the approximated discrete biharmonic distance<img src="./svgs/ecb2edc691a4d089460ad80188b5da8f.svg?invert_in_darkmode" align=middle width=13.9955178pt height=29.190975pt/>, exposing it through an Eigen-based API. The function

```C++
igl::biharmonic_distance_approx(V,F,k,D);
```

computes the approximated discrete biharmonic distances between each pair of vertices in V, using the first <img src="./svgs/63bb9849783d01d91403bc9a5fea12a2.svg?invert_in_darkmode" align=middle width=9.07536795pt height=22.8310566pt/> eigenvectors.



TODO: add image

Distance comparison between the exact approach (left) and the approximate approach (right).



TODO: update table

|                            | Exact | Approximate (with <img src="./svgs/071678d5a1ade4879c03e94c329c0314.svg?invert_in_darkmode" align=middle width=61.71225885pt height=22.4657235pt/>) |
| -------------------------- | ----- | ------------------------------------------------------------ |
| cactus (with ??? vertices) | ???   | ???                                                          |
| lucy (with ??? vertices)   | ???   | ???                                                          |

Runtime comparison between the exact approach and the approximate approach.



### References

1. Y. Lipman, R. M. Rustamov, T. A. Funkhouser, [Biharmonic Distance](https://www.cs.princeton.edu/~funk/biharmonic.pdf), 2010.

