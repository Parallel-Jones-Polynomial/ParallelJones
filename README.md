# ParallelJones
python package to evaluate Jones polynomial of curves in 3 space

© 2024 Kasturi Barkataki, Eleni Panagiotou

If you use this code, please cite the following paper:

**Kasturi Barkataki and Eleni Panagiotou. A parallel algorithm for the exact computation of the Jones polynomial. (In final stages of Preparation), 2025.**

Visit https://www.elenipanagiotou.com/ for updated information.

---

## Table of Contents
- [What is this package?](#what-is-this-package)
- [Getting Started](#getting-started)
- [Usage](#usage)
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
- **Anaconda Python interpreter Python 3.7+**
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

The  **ParJones.py** calculates the Jones polynomial of open or closed curves in 3D space.  It can be run from the command-line as follows :

```./anaconda3/bin/python ParJones.py num_projections input_knot closed parallel```
 

### **Input Arguments**  
| **Parameter**     | **Description**                              | **Type**
|--------------------|----------------------------------------------|-------
| `num_projections` | Number of projections for computation        | integer
| `input_knot`      | A numpy array of open/closed curves where each curve is given as a list of 3D coordinates.| string
| `closed`          | Curve type (`1` for closed, `0` for open)    | integer
| `parallel`        | Execution mode (`1` for parallel, `0` serial)| interger

---

### **Example Execution**  

The Jones polynomial of a closed trefoil in 3-space using the parallel algorithm over 1 projection:

./anaconda3/bin/python ParJones.py 1 '[[[1, 0, 0],[4, 0, 0],[1, 6, 2],[0, 2, -5],[5, 2, 5],[4, 6, -2]]]' 1 1


#### **Interpreting the Output**

The output of  is a dictionary where the **keys** represent the powers of the variable $A$ in the Jones polynomial, and the **values** correspond to the coefficients for each respective power of $A$.

```python
{-9.0: -0.0, -8.0: -0.0, -4.0: 1.0, -12.0: 1.0, -16.0: -1.0, 0: 0.0}
```
The above output implies that the Jones polynomial of the trefoil is $A^{-4}+A^{-12}-A^{-16}$.


## License
BSD 3-Clause "New" or "Revised" License
Copyright (c) 2021, Kasturi Barkataki and Eleni Panagiotou
All rights reserved.

---

## Contributors
Kasturi Barkataki and Eleni Panagiotou



