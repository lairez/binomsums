# This is the problem 11914 the American Mathematical Monthly (May 2016)
# The following sum is zero.

S := Sum(Sum(x^n*y^m*Sum((-4)^(-k)*Binomial(n-k,k-1)*Sum((-2)^(-j)*Binomial(n+1-2*k,j-1)*Binomial(m-k,3*m-j),j=1..3*m), k=1..n), n=1..infinity), m=1..infinity);

R, ord := BinomSums[sumtores](S, u);
#  -> 0, []


