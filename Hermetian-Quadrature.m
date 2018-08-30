f[x_] = Sin[x] E^-x;
f'[x_] = D[f[x], x];
partition = Range[0, Pi, (Pi/4)];

polynomial[int_, var_] := 
  With[{i = int, x = var}, 
   prod = Product[(x - partition[[j]])/(partition[[i]] - 
        partition[[j]]), {j, Complement[Range[1, 5], {i}]}]; 
   Return[prod];];

H[int_, var_] := With[{i = int, x = var},
   p'[i_, x_] = D[polynomial[i, x], x];
   (polynomial[i, x]^2) (1 - 
      2*p'[i, partition[[i]]]*(x - partition[[i]]))];

S[var_, x_] := 
  With[{i = var}, (polynomial[i, x])^2 (x - partition[[i]])];

Print["Hermetian Interpolant"];
hermiteInterpolation = 
 Sum[f[partition[[i]]]*H[i, x], {i, 1, 5}] + 
  Sum[f'[partition[[i]]]*S[i, x], {i, 1, 5}]

Print["Hermite Quadrature"]
hermiteQuadrature = 
 Sum[f[partition[[i]]]*Integrate[H[i, x], {x, 0, Pi}], {i, 1, 5}] + 
   Sum[f'[partition[[i]]]*Integrate[S[i, x], {x, 0, Pi}], {i, 1, 
     5}] // N

Print["NIntegrate"]
NIntegrate[f[x], {x, 0, Pi}]