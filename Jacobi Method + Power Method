n = 50;
a = DiagonalMatrix[Array[.25 &, n - 1], -1];
d = DiagonalMatrix[Array[.5 &, n]];
c = DiagonalMatrix[Array[.25 &, n - 1], 1];
matrixA = a + d + c;
b = 1/Array[# &, n];

(*#1*)Print["Jacobi method"];
abar = matrixA - d;
x = Array[1. &, n];
list1 = {};
i = 1;
For[i > 1, Norm[ b - matrixA.x] > 10^-5, i++,
 var = (Inverse[d].(-abar)).x + Inverse[d].b;
 x = var;
 ]; Jacobi = var;
AppendTo[list1, {"Final Iteration: " <> ToString[i], 
   " Residual value: ", Norm[b - matrixA.x], " Solution: ", var}];
list1

(* #2 Power Method *)
n = 50;
a = DiagonalMatrix[Array[.25 &, n - 1], -1];
d = DiagonalMatrix[Array[.5 &, n]];
c = DiagonalMatrix[Array[.25 &, n - 1], 1];
matrixA = a + d + c;
b = 1/Array[# &, n];
x0 = d.b;
var = b - matrixA.x0;
power = x0;
list = {};
j = 1;
For[i = 1, i < 1000 && j > 10^-5, i++,
  e = b - matrixA.(x0 + d.var);
  j = (Norm[e])/(Norm[b]);
  power = (power.matrixA)/Norm[power.matrixA];
  ];
AppendTo[list, "Spectral Radius: " <> ToString[Max[power]]];
list