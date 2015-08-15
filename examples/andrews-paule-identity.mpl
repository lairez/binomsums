# Andrews-Paule identity. See §7.1 in “Multiple binomial sums”.

S := Sum(z^n*Sum(Sum(Multinomial([i,j])^2*Multinomial([2*n-2*j,2*n-2*i]), i=0..infinity), j=0..infinity), n=0..infinity);

R, ord := BinomSums[sumtores](S, u);

# Checking
BinomSums[rser](R, ord, 6);
BinomSums[computesum](S, 5);

# We obtain the following differential operator annihilating the generating function S:
dop := 1048576*D^6*z^8 + 2883584*D^6*z^7-40960*D^6*z^6 +
22020096*D^5*z^7-29696*D^6*z^5 + 64225280*D^5*z^6 + 1296*D^6*z^4 +
1269760*D^5*z^5 + 146407424*D^4*z^6-605184*D^5*z^4 + 455041024*D^4*z^5 +
16848*D^5*z^3 + 24412672*D^4*z^4 + 363069440*D^3*z^5-3632352*D^4*z^3 +
1211465728*D^3*z^4 + 59292*D^4*z^2 + 106845184*D^3*z^3 +
305827840*D^2*z^4-7352832*D^3*z^2 + 1109626112*D^2*z^3 + 58320*D^3*z +
139138736*D^2*z^2 + 60272640*D*z^3-4247073*D^2*z + 244005120*D*z^2 + 9720*D^2 +
42117840*D*z + 691200*z^2-374625*D + 3369600*z + 996300;

# Conversion to a differential equation with initial condition
init := BinomSums[computesum](S, 6):
deq := {
  DETools[diffop2de](dop, y(z), [D,z]),
  seq( (D@@i)(y)(0) = i!*coeff(init, z, i), i=0..4 )
  };

# Conversion to a recurrence.
rec := gfun[diffeqtorec](deq, y(z), u(n));

# Resolution of the recurrence.
rsolve(rec, u(n));

