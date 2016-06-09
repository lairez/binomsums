# Strehl's identity. See §7.2.1 in “Multiple binomial sums”.

# The identity is S1 = S2.
S1 := Sum(z^n*Sum(Binomial3(n,k)*Multinomial([n,k])*Sum(Binomial3(k,j)^3, j=0..infinity), k=0..infinity), n=0..infinity);
S2 := Sum(z^n*Sum(Binomial3(n,k)^2*Binomial3(n+k,k)^2, k=0..infinity), n=0..infinity);

# Quick check
BinomSums[computesum](S1, 5) - BinomSums[computesum](S2, 5);

R1, ord1 := BinomSums[sumtores](S1, u);

# Consistency check
BinomSums[rser](R1, ord1, 6);

R2, ord2 := BinomSums[sumtores](S2, u);

# Consistency check
BinomSums[rser](R2, ord2, 6);

# For R1 and R2 we find the same Picard Fuchs equation :
# D^3 + (6*t^2 - 153*t + 3)/(t^3 - 34*t^2 + t)*D^2 + (7*t^2 - 112*t + 1)/(t^4 - 34*t^3 + t^2)*D + (t - 5)/(t^4 - 34*t^3 + t^2)
# Since the initial conditions are the same, S1 = S2.
