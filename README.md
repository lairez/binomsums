BinomSums
=========

*BinomSums* is a package, for the computer algebra system Maple 17 or later,
providing functions to handle multiple binomial sums.  It implements the
algorithms described in the paper (soon on the arXiv):

> A. Bostan, P. Lairez, and B. Salvy, “Multiple binomial sums”

Table of Contents
-----------------

* [Warning](#warning)
* [License](#license)
* [Installation](#installation)
* [Usage](#usage)
  * [Binomial sums](#binomial-sums)
  * [Examples of binomial sums](#examples-of-binomial-sums)
  * [Description of the main functions](#description-of-the-main-functions)
    * [sumtores](#sumtores)
    * [ratres](#ratres)
    * [geomred](#geomred)
    * [computesum](#computesum)
    * [rser](#rser)
  * [Other functions](#other-functions)
    * [geomredall](#geomredall)
    * [sumtoct](#sumtoct)
    * [addnewgf](#addnewgf)


Warning
-------

*This is a very preliminary implementation.*
It is likely to contain bugs. It is provided *as is*.
Do not hesitate to submit bug reports of pull requests.


License
-------

*BinomSums* is released under the terms of the MIT license.

See the file LICENSE for more information.


Installation
------------

1. Execute `make`, this produces the file `binomsums.mla`.
2. Check that the variable `libname` in Maple contains the path to `binomsums.mla` or to its parent directory.
3. Load the package in Maple with `with(BinomSums);`


Usage
-----

### Binomial sums

A *binomial sum* is a Maple expression made with the following rules.  For
simplicity, let us assume that Maples indeterminates split in two sets: the
continuous variables and the discrete ones. An `<integer form>` is a polynomial
of degree one in the discrete variables, with integer coefficients. 

A *supergeometric* sequence `<supergeom>` can be:

* A polynomial in the discrete variables
* `<rational function of continous variables>^<integer form>`
* `<supergeom> + <supergeom>`
* `<supergeom> * <supergeom>`

A `<binomial sum>` can be either *primitive*:

* `Binomial(<integer form >, <integer form>)`
* `Delta(<integer form >)`
* `H(<integer form >)`
* `Catalan(<integer form >)`
* `Multinomial([<integer form >, ...])`
* `<supergeom>`

Or a `<binomial sum>` can be *composed*:

* `<binomial sum> + <binomial sum>`
* `<binomial sum> * <binomial sum>`
* `<binomial sum>^<integer>`
* `Sum(<binomial sum>, <discrete variable>=<integer form>..(<integer form> || infinity))`

A generating function `<gfun>` is a binomial sum whose discrete variables are
all bound by a `Sum`.

### Examples of binomial sums

```
S1 := Binomial(n,k)^2*Multinomial([n,k]);

S2 := Sum(t^n*Sum(Binomial(n,k)^2*Multinomial([n,k]), k=0..infinity), n=0..infinity);

S3 := Sum(Sum(x^n*y^m*Binomial(n+2*m, m), n=0..infinity), m=0..infinity);

S4 := Sum(t^n*Sum(Binomial(n,n-k), k=0..infinity), n=0..infinity);
```

Note that `S2`, `S3` and `S4` are generating functions.

### Description of the main functions

#### `sumtores`

`sumtores(S :: <gfun>, v :: name, [geomred = false]) -> ratpoly, list(name)`

Returns a rational function `R` and a list of names `ord`.

The variables in `R` are the continuous variables in `S` and `v[1]`, `v[2]`,
etc.

The generating function `S` is the residue of `R` with respect to the `v[i]`'s.
The list `ord` gives, by ascending order, the order of the variables defining
the iterated Laurent series field in which the residue is defined (see the
paper for definition).

`sumtores` fails if it does not find an order on the variables that makes the
infinite sums converge.

If `geomred=false` is given, then the *geometric reduction* (see the paper for
definition) is not performed.

##### Examples

```
> sumtores(S2, v);
                                 1
  ---------------------------------------------------------------, [t, v[1], v[2]]
  (v[1] v[2] - v[1] - v[2]) t + v[1] v[2] (-1 + v[2]) (-1 + v[1])

> sumtores(S3, v);
                                v[1]
                  ---------------------------------, [x, y, v[1]]
                  (y + v[1] (-1 + v[1])) (x - v[1])

> sumtores(S4, v);
Error, (in solvecons) inconsistent

```

#### `ratres`

```
ratres(R :: ratpol, v :: name, ord :: list(name)) -> ratpoly or FAIL
```

Returns `FAIL` or a rational function `T` that is the residue of `R` with
respect to the variable `v` in the iterated Laurent series field defined by the
ordering `ord`. It implements Algorithm 3 of the paper.

##### Example

```
> ratres(1/(t*(1+z)^19+z-t), z, [t, z]);
                                      1
                                   --------
                                   19 t + 1
```

#### `geomred`

```geomred(R :: ratpoly, ord :: list(name), params :: set(name)) -> ratpoly```

Tries to apply `ratres` with all the variables of `ord` that are not in
`params`, one after the other, and returns the result.

So it returns a rational function `T` with less variables, hopefully, than `R`
such that its residue with respect to the variables in `R` that are not in
`params` is equal to the residue of `S` with respect to the variables in `S`
that are not in `params`.

#### `computesum`

`computesum(S :: <binomial sum>, maxn :: integer)`

Replace `infinity` by `maxn` in `S`, and then compute the sum.  Useful for
checking that things are consistent.

##### Example
```
> computesum(S2, 6); 
                   6          5         4        3       2
           104959 t  + 11253 t  + 1251 t  + 147 t  + 19 t  + 3 t + 1
```

#### `rser`

`rser(R :: ratpoly, vars :: list(name), n :: posint) -> truncated power series`

Compute the first `n` terms of the power series expansion of the residue of `R`
with respect to all the variables of `vars` except the first one.  The ordering
of `vars` defines the iterated Laurent series field in which the computation
takes place.  Useful for checking that things are consistent.

##### Example
```
> rser(sumtores(S2, v), 7);
                         2        3         4          5           6
           1 + 3 t + 19 t  + 147 t  + 1251 t  + 11253 t  + 104959 t
```


### Other functions

#### `geomredall`

`geomredall(R :: ratpoly, ord :: list(name), params :: set(name)) -> set(ratpoly)`

The same as `geomred` but tries every possible order to eliminate the
variables.

#### `sumtoct`

`sumtoct(S :: <binomial sum>, v :: name) -> <supergeom>` 

Returns `T :: <supergeom>`.

If `S` contains no infinite sum, the continuous variables of `T` are those of
`S` and `v[1]`, `v[2]`, etc.  In this case, for any integer value of the
discrete variables, the value of `S` is the constant term of `T`, with respect
to the variables `v[i]`'s, in any order.

`T` will contain an extra variable `_W[_k]` for each infinite sum
`k=a..infinity`. They are used internally to determine an order on the
variables that makes things converge.

##### Examples

```
> sumtoct(S1, v);                                                                   
             (k + 1)     (n - k)           (k + 1)     (n - k)
1/((1 - v[1])        v[1]        (1 - v[2])        v[2]        (1 - v[3] - v[4])

        n     k
    v[3]  v[4] )

> sumtoct(S2, v);
- v[4] v[1] v[2] v[3]/((-v[1] v[2] v[3] + t _W[_n])

    (_W[_k] v[1] v[2] - v[1] v[2] v[4] + v[1] v[4] + v[2] v[4] - v[4])

    (-1 + v[3] + v[4]))

```

#### `addnewgf`

*TODO*



