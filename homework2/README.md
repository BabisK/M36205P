50 runs:
Total time recursive: 2.323047
Total time LAPACK LU inverse: 3.207433

Once: Recursive first (cache invalidated)
Total time recursive: 0.078320
Total time LAPACK LU inverse: 0.082500

Once: LAPACK first (cache invalidated)
Total time recursive: 0.057625
Total time LAPACK LU inverse: 0.109023

R:
> i1
   user  system elapsed 
  0.043   0.000   0.044 
> i2
   user  system elapsed 
   0.76    0.00    0.76

gcc -O3 recinvert.c -lcblas -llapacke 
