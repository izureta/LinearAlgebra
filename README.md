# LinearAlgebra
This is a simple c++ linear algebra library

# Implemented
* Matrix class (with rational numbers by default)
* Rational numbers class
* Some linear algebra algorithms: Gauss algorithm, Jacobi algorithm (of diagonalizing symmetrical bilinear form), Gram-Schmidt process

# Examples
* Basic matrix operations

```cpp
    >>> #include "LinearAlgebra.h"

    >>> Matrix A({{1, 2, 3, 4}});
    >>> Matrix B({{-1}, {-2}, {-3}, {-4}});
    >>> (B * (A + A)).Print();
    -2 -4 -6 -8
    -4 -8 -12 -16
    -6 -12 -18 -24
    -8 -16 -24 -32
```

* Solving system of linear equations

<img src="https://latex.codecogs.com/svg.image?&space;\begin{cases}&space;&space;&space;x&plus;y=3&space;\\&space;&space;&space;x&plus;z=4&space;\\&space;&space;&space;y&plus;z=0&space;\end{cases}" />

```cpp
    >>> Matrix A({{1, 1, 0, 3},
            {1, 0, 1, 4},
            {0, 1, 1, 0}});
    >>> A.Gauss().Print();
    1 0 0 7/2
    0 1 0 -1/2
    0 0 1 1/2
```

<img src="https://latex.codecogs.com/svg.image?&space;\begin{cases}&space;&space;&space;x=\frac{7}{2}&space;\\&space;&space;&space;y=\frac{-1}{2}&space;\\&space;&space;&space;z=\frac{1}{2}&space;\end{cases}" />

* Quadratic form normalization

<img src="https://latex.codecogs.com/svg.image?f(x)&space;=&space;x_1^2&plus;2x_1x_2&plus;2x_2^2&plus;2x_2x_3&plus;x_3^2" />

Correspondant symmetrical matrix

<img src="https://latex.codecogs.com/svg.image?B=\begin{pmatrix}&space;1&space;&&space;1&space;&&space;0&space;\\&space;1&space;&&space;2&space;&&space;1&space;\\&space;0&space;&&space;1&space;&&space;1&space;\end{pmatrix}" />

```cpp
    >>> Matrix B({{1, 1, 0},
            {1, 2, 1},
            {0, 1, 1}});
    >>> B.Jacobi().first.Print();
    1 0 0
    0 1 0
    0 0 0
```

<img src="https://latex.codecogs.com/svg.image?f(x)=g(y)=y_1^2&plus;y_2^2" />

* Matrix with complex numbers

<img src="https://latex.codecogs.com/svg.image?B=\begin{pmatrix}&space;1&0&i&space;\\&space;0&i&0&space;\\&space;0&2&0&space;\end{pmatrix}" />

```cpp
    >>> #include <complex>
    >>> using cmpl = std::complex<int64_t>;
    >>> Matrix<cmpl> B({{cmpl(1, 0), cmpl(0, 0), cmpl(0, 1)},
                        {cmpl(0, 0), cmpl(0, 1), cmpl(0, 0)},
                        {cmpl(0, 0), cmpl(2, 0), cmpl(0, 0)}});
    >>> std::cout << B.Rank();
    2
```

<img src="https://latex.codecogs.com/svg.image?rank(B)&space;=&space;2" />
