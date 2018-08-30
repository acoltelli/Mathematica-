(*#1 *)

Print["Gauss Seidel"];
lower = d + 
   Table[If[(i > j), matrixA[[i, j]], 0], {i, 1, n}, {j, 1, n}];
upper = Table[If[i < j, matrixA[[i, j]], 0], {i, 1, n}, {j, 1, n}];
x = Array[1. &, n];
list2 = {};
i = 1;
For[i > 1, Norm[ b - matrixA.x] > 10^-5, i++,
 var = (Inverse[lower].-upper).x + Inverse[lower].b;
 x = var;]; Gauss = var;
AppendTo[list2, {"Final Iteration: " <> ToString[i], 
   " Residual value:" , Norm[ b - matrixA.x], 
   " Solution: " <> ToString[var]}];
list2




(*#2 *)
n = 100;
a = DiagonalMatrix[Array[.25 &, n - 1], -1];
d = DiagonalMatrix[Array[.5 &, n]];
c = DiagonalMatrix[Array[.25 &, n - 1], 1];
matrixA = a + d + c;
b = 1/Array[# &, n];
list = {};
Print["Gauss Seidel(100x100 matrix)"];
lower = d + 
   Table[If[(i > j), matrixA[[i, j]], 0], {i, 1, n}, {j, 1, n}];
upper = Table[If[i < j, matrixA[[i, j]], 0], {i, 1, n}, {j, 1, n}];
x = Array[1. &, n];
i = 1;
For[i > 1, Norm[ b - matrixA.x] > 10^-5, i++,
 var = (Inverse[lower].-upper).x + Inverse[lower].b;
 x = var;]; Gauss = var;
AppendTo[list, {"Final Iteration: " <> ToString[i], 
   " Residual value:" , Norm[ b - matrixA.x], 
   " Solution: " <> ToString[var]}];
list