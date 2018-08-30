n = 100;
alpha = .5;
dx = 0.1;
dt = 0.01;
lambda = (alpha*dt)/dx^2;

matrixA = 
  DiagonalMatrix[Array[-lambda &, n - 1], -1] + 
   DiagonalMatrix[Array[((2 + 2 lambda)) &, n]] + 
   DiagonalMatrix[Array[-lambda &, n - 1], 1];
identity = IdentityMatrix[n];
matrixA[[1]] = identity[[1]];
matrixA[[n]] = identity[[n]];
matrixA = .5*matrixA;

matrixB = 
  DiagonalMatrix[Array[lambda &, n - 1], -1] + 
   DiagonalMatrix[Array[((2 - 2 lambda)) &, n]] + 
   DiagonalMatrix[Array[lambda &, n - 1], 1];
identity = IdentityMatrix[n];
matrixB[[1]] = identity[[1]];
matrixB[[n]] = identity[[n]];
matrixB = .5*matrixB;

u1 = Array[20 &, n];
u1[[1]] = 0;
u1[[n]] = 0;

states = {};
AppendTo[states, u1];
Do[u = Inverse[matrixA].matrixB.states[[i]]; 
  AppendTo[states, u], {i, 1000}];

p1 = ListLinePlot[states[[100]], PlotLabels -> Placed[{"t=10"}, Above]]; 
p2 = ListLinePlot[states[[200]], PlotLabels -> Placed[{"t=9"}, Above]]; 
p3 = ListLinePlot[states[[300]], PlotLabels -> Placed[{"t=8"}, Above]]; 
p4 = ListLinePlot[states[[400]], PlotLabels -> Placed[{"t=7"}, Above]]; 
p5 = ListLinePlot[states[[500]], PlotLabels -> Placed[{"t=6"}, Above]];
p6 = ListLinePlot[states[[600]], PlotLabels -> Placed[{"t=5"}, Above]]; 
p7 = ListLinePlot[states[[700]], PlotLabels -> Placed[{"t=4"}, Above]]; 
p8 = ListLinePlot[states[[800]], PlotLabels -> Placed[{"t=3"}, Above]];
p9 = ListLinePlot[states[[900]], PlotLabels -> Placed[{"t=2"}, Above]];
p10 = ListLinePlot[states[[1000]], PlotLabels -> Placed[{"t=1"}, Above]];
Show[p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, 
 PlotLabel -> "Crank Nicolson"]