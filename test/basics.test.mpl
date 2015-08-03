with(BinomSums):
with(CodeTools):

### rser
Test(
  coeff(rser(1/(1-2*t), [t], 10),t,8),
  2^8);

Test(
  coeff(rser(1/(1-2*t)/u[1], [t,u[1]], 10),t,8),
  2^8);

Test(
  coeff(rser(1/(1-2*t)/u[1]/(u[2]-u[1]), [t,u[1],u[2]], 10),t,8),
  2^8);

### packvars
Test(
  BinomSums[packvars]([ 1/(x[4]+x[6]), x[7] ], x),
  [{x[4] = x[1], x[6] = x[2]}, {x[7] = x[3]}]
);


### geomsum
Test(
  BinomSums[geomsum](1, k=0..n),
  n+1
);

Test(
  BinomSums[geomsum](k, k=0..n),
  n*(n+1)/2,
  `expand`
);

# when the sum is infinite, geomsum add a variable
Test(
  BinomSums[geomsum](k, k=0..infinity),
  _W[_k]/(1-_W[_k])^2,
  `normal`
);

# when the sum is infinite on the left, geomsum yields an error
Test(
  BinomSums[geomsum](k, k=-infinity..0),
  "",
  `testerror`
);

# it should recognize when the input does not really depend on the summation variable
Test(
  BinomSums[geomsum]((-1+u[1])^(-n-1+j)*(-1+u[2])^(-n+j-1)*(1-u[1])^(-j)*(1-u[2])^(-j)*u[1]^(-n)*u[2]^(-n), j=0..n),
  u[1]^(-n)*u[2]^(-n)*(-1+u[1])^(-n-1)*(-1+u[2])^(-n-1)*(n+1),
  `normal`
);

