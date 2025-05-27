# ParallelJones
python package to evaluate Jones polynomial of curves in 3 space

© 2024 Kasturi Barkataki, Eleni Panagiotou

If you use this code, please cite the following paper:

**Kasturi Barkataki and Eleni Panagiotou. A parallel algorithm for the exact computation of the Jones polynomial. (In final stages of Preparation), 2024.**

Visit https://www.elenipanagiotou.com/ for updated information.

---

## Table of Contents
- [What is this package?](#what-is-this-package)
- [Getting Started](#getting-started)
- [Usage](#usage)
- [Overview of Jones.py](#overview-j)
- [Overview of functions.py](#overview-f)
- [License](#license)
- [Contributors](#contributors)

---

## What is this package?

This package computes the Jones Polynomial of one or more open/closed curves (knots or links) in 3-space. The package supports serial and parallel execution.
The input can be specified as 3D coordinates of the curves, and the output is the Jones Polynomial, which is a topological invariant for closed knots/links and
a continuous function of curve coordinates for collections of open curves.

---

## Getting Started <a name="getting-started"></a>

To use this package, the following dependencies should be installed:

### Dependencies:
- **Python 3.7+**
- **Required Python packages**:
  - `numpy`: For numerical computations and array manipulations.
  - `math`: Provides mathematical functions (e.g., trigonometric, logarithmic).
  - `copy`: For creating deep copies of data structures.
  - `multiprocessing`: To enable parallel computation of the Jones Polynomial.
  - `time`: For measuring the execution time.
  - `itertools`  : For operations like creating combinations, permutations
  - `sys`: For memory management, such as checking the size of objects.
  - **functions.py**: A custom module that contains essential functions required for calculating the Jones Polynomial. 
---
## **Usage**

The  **Jones.py** calculates the Jones polynomial of open or closed curves in 3D space. It takes the following inputs.  

### **Inputs**  
| **Parameter**     | **Description**                              |
|--------------------|----------------------------------------------|
| `num_projections` | Number of projections for computation        | 
| `input_knot`      | A numpy array of open/closed curves where each curve is given as a list of 3D coordinates.| 
| `closed`          | Curve type (`1` for closed, `0` for open)    |
| `parallel`        | Execution mode (`1` for parallel, `0` serial)| 

---

### **Example Execution**  

```python
import numpy as np  
from functions import *
from Jones import *

# User-defined inputs
num_projections = 1  
trefoil = np.array([[[1, 0, 0], [4, 0, 0], [1, 6, 2], [0, 2, -5], [5, 2, 5], [4, 6, -2]]], dtype=object)  
closed = 1       # Closed curve  
parallel = 1     # Parallel execution  

# Run Jones Polynomial calculation
jones_execution(num_projections, trefoil)
```


#### **Interpreting the Output**

The output of `jones_execution` is a dictionary where the **keys** represent the powers of the variable $A$ in the Jones polynomial, and the **values** correspond to the coefficients for each respective power of $A$.

```python
{-9.0: -0.0, -8.0: -0.0, -4.0: 1.0, -12.0: 1.0, -16.0: -1.0, 0: 0.0}
```
The above output implies that the Jones polynomial of the trefoil is $A^{-4}+A^{-12}-A^{-16}$.


---
## Overview of Jones.py

The user calls the `jones_execution` function to compute the Jones polynomial. This function is built by means of a collection of
internal functions which are briefly described in the following table.

### Internal Functions
| **Function**                          | **Purpose**                                                                                                    |
|---------------------------------------|----------------------------------------------------------------------------------------------------------------|
| `get_jones_poly(coords, proj_vector, Inds)`  | Computes the Jones polynomial by projecting 3D coordinates onto a 2D plane and calculating the polynomial.   |
| `crossing_matrices_bool(proj, inds, depth)` | Generates crossing matrices and boolean masks representing edge crossings, over/under status, and orientations. |
| `get_partial_poly(bool_mask, over_or_under, right_or_left, inds, clos_perm, modf_inds, states, S)` | Recursively computes the bracket polynomial using crossing matrices and edge relations.                          |
| `split_jones(arg_1, arg_2, clos_perm, Wr)`    | Splits the knot/link into smaller linkoids for parallel computation of the Jones polynomial.                  |
| `find_arc(bool_mask, inds)`                | Finds arcs in the crossing matrices to split the knot or link into linkoids.                                  |
| `get_writhe(proj, inds)`                    | Calculates the writhe (topological invariant) of the knot or link.                                            |
| `form_indices(coords, proj_vector)`        | Forms indices that define the knot or link structure based on the 3D coordinates and the projection vector.    |
| `process_projection(coords, proj_vector)`   | Projects the 3D coordinates onto a 2D plane based on the specified projection vector.                         |


---
## Overview of functions.py
Utility functions for **knot theory** computations used in **Jones.py**, including Reidemeister moves, graph simplification, and polynomial operations.

| **Function**               | **Description**                                                   |
|----------------------------|-------------------------------------------------------------------|
| `get_random_proj()`   | Generates a random point on the surface of a unit sphere to specify a direction of projection.        |
| `fibonacci_sphere(samples)` | Uniformly samples points on a sphere using the Fibonacci lattice. |
|`get_two_vec(proj_vec)`| Finds two orthonormal vectors orthogonal to a given 3D vector.   |
| `getRoots(aNeigh)`    | Finds connected components in a graph.                           |
| `Loops(inds, clos_perm)` | Partitions vertices into disjoint sets based on closure relations. |
| `det2x2(A)`            | Computes the determinant of a 2x2 matrix.                        |
|`Cramer(A)`           | Solves a 2x2 linear system using Cramer’s rule.                  |
|`find_conn(d, g, inds, bool_mask)` | Checks if two vertices in a graph are connected.           |
|`simplification(BM, over_or_under, inds)` | Simplifies a knot graph using Reidemeister moves.        |
|`max_len(Ch)`          | Finds the maximum length of elements in a list.                 |
|`J_mult(P, Q)`        | Multiplies two polynomials.                                      |
|`J_add(P, Q)`         | Adds two polynomials.                                            |
|`J_smult(a, P)`        | Scales a polynomial by a scalar.                                |
|`dfactor(N)`          | Returns the polynomial $(-A^{-2}-A^{2})^{N}$.             |

---

## License
BSD 3-Clause "New" or "Revised" License
Copyright (c) 2021, Kasturi Barkataki and Eleni Panagiotou
All rights reserved.

---

## Contributors
Kasturi Barkataki and Eleni Panagiotou



