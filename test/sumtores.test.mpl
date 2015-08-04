with(BinomSums):
with(CodeTools):

### sumtoct

S := Sum(Binomial(n,k),k=0..n);
R := sumtoct(S, u);

Test(residue(subs(n=10,R)/u[1], u[1]=0), 2^10);


### sumtores

S := Sum(t^n*Sum(Binomial(n,k),k=0..n), n=0..infinity);
R, ord := sumtores(S, u);

Test(nops(ord), 1);
Test(R, 1/(1-2*t), `normal`);

