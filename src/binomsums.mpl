BinomSums := module()

option package;

local
  gfdict, gfdictcomp, gfnames, multinomial, packvars, geomsum,
  sumtores1, solvecons, linearorder, slopes, asylt, inorout_base, ratres2,
  geomsum_, simpfacts, simpfacts0 ;

export
  sumtores, SimpleImpl, addnewgf, isconvergent, rser, computesum, asy, inorout,
  hermitered, ratres, geomred, geomredall, sumtoct;

  # This module provides a simple implementation of the computation of integral
  # representation of binomial sums.
  SimpleImpl := module()
  
    export sumtores;

    # Input :
    #   - U a binomial sum, as a Maple expression
    #   - v a name
    #   - num an integer (default: 1)
    # Output :
    #   A rational function R(v[num],...,v[num+k]) such that U = res R.
    sumtores := proc(U, v :: name, num :: integer := 1)
      local L, first, rest, rat_first, num_first, i;
      if type(U, `+`) then
        return normal(map(sumtores1, U, v, num));
      elif type(U, `*`) or type(U, `^`(anything, posint)) then
        first := op(1, U);   # the first factor
        if type(U, `*`) then
          rest := subsop(1=1, U);  # the rest of the product
        else
          rest := first^(op(2, U) - 1);
        end if;
        rat_first := sumtores1(first, v, num);
        num_first := max(num, op(map(op, indets(rat_first, specindex(v)))));
        return normal(rat_first*sumtores1(rest, v, num_first + 1));
      elif type(U, specfunc(Delta)) then
        return v[num]^op(U);
      elif type(U, specfunc(Binomial)) then
        return (1+v[num])^op(1, U)/v[num]^(op(2, U)+1); 
      elif type(U, specfunc(Sum)) then
        return normal(sum(expand(sumtores1(op(1, U), v, num)), op(2, U)));
      else
        return U/v[num];
      end if;
    end proc;

  end module:



#### SUM TO PERIOD

gfdict := table();
gfdictcomp := table();
gfnames := {};

# Input :
#   - name, a name
#   - a function gf such that name(n) stands for res( gf(n) )
#   - a function comp such that comp( n ) = res( gf(n) )
#
# Effect :
#   Add 'name' as a basic block for constructing binomial sums
addnewgf := proc(name :: name, gf, comp)
  #global gfdict, gfdictcomp, gfnames;
  gfdict[name] := gf;
  gfdictcomp[name] := comp;
  gfnames := gfnames union {name};
end proc;

addnewgf(Binomial,
  ((v,n,k) -> 1/(1-v[1])^(k+1)/v[1]^(n-k)),
  binomial);

addnewgf(Binomial2,
  ((v,n,k) -> (1+v[1])^n/v[1]^(k)),
  binomial);

addnewgf(Multinomial,
  ((v,L) -> 1/(1-add(v[i],i=1..nops(L)))/mul(v[i]^(L[i]),i=1..nops(L))),
  (L ->multinomial(L)));

addnewgf(Catalan,
  ((v,n) -> (1+v[1])^(2*n)*(1-v[1])/v[1]^(n)),
  (n -> 1/(n+1)*binomial(2*n,n) ));

addnewgf(H,
  ((v,n) -> v[1]^(-n)/(1-v[1])),
  n -> `if`(n >= 0, 1, 0));

addnewgf(Delta,
  ((v,n) -> v[1]^n),
  n -> `if`(n = 0, 1, 0));

#addnewgf(CT,
#  ( (R, v) -> CT(R,v) ),
#  (R, v) -> residue(R/v,v=0));

multinomial := proc( x :: list(numeric))
  if nops(x) <= 1 then
    return 1;
  elif x[1] < 0 then
    return 0;
  else
    return binomial(`+`(op(x)), x[1])*multinomial(x[2..-1]);
  end if;
end;


# Outputs a list of sets s[1],...,s[nops(L)] such that the substitutions
# subs(s[i], L[i]) guarantees that the variables of the elements of L does not
# overlap and and are numbered consecutively.
#
# Example : packvars([ 1/(x[4]+x[6]), x[7] ], x);
#    => [{x[4] = x[1], x[6] = x[2]}, {x[7] = x[3]}]
packvars := proc(L :: list, v :: name)
  local Lv, Ls, i, j, l;
  Lv := map2(select, has, map(indets, L, name), v);
  Ls := NULL;
  i := 0;
  for l in Lv do
    Ls := Ls, {seq(l[j]=v[i+j], j=1..nops(l))};
    i := i + nops(l);
  end do;
  return [Ls];
end proc;

#
# (Alternative implementation)
geomsum_ := proc(S, bounds)
  local T, svar, infb, supb, ret, prim:
  T := expand(normal(S));
  if type(T, `+`) then
    return normal(map(geomsum_, T, bounds));
  fi;
  
  svar := op(1, bounds);
  infb := op([2, 1], bounds);
  supb := op([2, 2], bounds);

  ASSERT(infb <> -infinity);
  if supb = infinity then
    T := _W[_||svar]^svar*T;
  fi;

  prim := SumTools[Hypergeometric][Gosper](T, svar);
  ret := -subs(svar=infb, prim);
  
  if supb <> infinity then 
    ret := ret + subs(svar=supb+1, prim);
  fi;
  
  return normal(ret);
end proc;

# Input :
#   - S, an expression of the form P(k)*A^k, where P is a polynomial
#   - bounds, an expression in the form k=a..b, where b can be infinity
#
# Returns an expression T without k such that
#   T = sum(S(k), k=a..b)
#
# If the upper bound b is `infinity' then it marks the sum
# with the extra variable _W[_k].
geomsum := proc(S, bounds)
  local T, svar, infb, supb, ret, prim, den, opT, denother, dendep:
  
  svar := op(1, bounds);
  infb := op([2, 1], bounds);
  supb := op([2, 2], bounds);

  T := normal(S);
  
  if infb = -infinity or supb = -infinity then
    error "Summation bounds cannot be -infinity.";
  fi;

  if supb = infinity then
    T := _W[_||svar]^svar*T;
  fi;

  # dendep (resp. denother) contains the factors of the denominator of T that
  # depends (resp. do not depend) on the summation variable.
  den := factor(denom(T));
  if type(den, `*`) then
    dendep, denother := selectremove(has, den, svar);
  elif has(den, svar) then
    dendep, denother := den, 1;
  else
    dendep, denother := 1, den;
  end if;

  T := expand(numer(T));
  if type(T, `+`) then
    opT := [op(T)];
  else
    opT := [T];
  end if;
  opT := map(`*`, opT, 1/dendep);
  
  # Factor simplification, to rewrite things like
  # A := -(-1+u[1])^(-n-1+j)*(-1+u[2])^(-n+j-1)*(1-u[1])^(-j)*(1-u[2])^(-j)*u[1]^(-n)*u[2]^(-n)
  # which actually does not depend on j.
  # On such an input, SumTools[IndefiniteSummation](A, j) does not work.
  opT := map(simpfacts, opT);
  opT := map(SumTools[IndefiniteSummation], opT, svar);

  ret := -convert(subs(svar=infb, opT),`+`);
  
  if supb <> infinity then 
    ret := ret + convert(subs(svar=supb+1, opT), `+`);
  fi;
  
  return (ret/denother);
end proc;


# Input:
#   - P, a product of expressions of the form F^a, where `F' is a polynomial
#   and `a' an expression that may depend on variables.
#
# Returns an equivalent product (when the variables in the exponent are
# integers) with the garantee that for any two factors F^a and G^b, F is not
# proportional to G.
simpfacts := proc(P)
  local t, i, ret;
  t := simpfacts0(P);
  ret := 1;
  for i in indices(t) do
    if op(i) = -1 then
      ret := ret*op(i)^(t[op(i)] mod 2);
    else
      ret := ret*op(i)^t[op(i)];
    end if;
  end do;
end proc;

# Input:
#   - same input as simpfacts
#
# Return an associative array T such that
# P = mul(F^a, (F, a) in T)
simpfacts0 := proc(P)
  local t, t1, ind, f, i, lc;
  if type(P, `^`) then
    t := simpfacts0(op(1, P));
    for f in indices(t) do
      t[op(f)] := t[op(f)]*op(2, P);
    end do;
    return t;
  elif type(P, `*`) then
    t := table();
    ind := {};
    for f in [op(P)] do
      t1 := simpfacts0(f);
      for i in indices(t1) do
        if not (i in ind) then
          ind := ind union {i};
          t[op(i)] := 0;
        end if;
        t[op(i)] := t[op(i)] + t1[op(i)];
      end do;
    end do;
    return t;
  else
    f := normal(expand(P));
    lc := lcoeff(f);
    f := normal(expand(f/lc));
    return table([ lc = 1, f = 1 ]);
  end if;
end proc;

# Given a binomial sum, returns an expression such that the binomial sum equal
# the constant term of the expression.
#
# The output may contains extravariables _W[_k] to track infinite summations
sumtoct := proc(S, v :: name)
  local L;
  #global gfdict, gfnames;
  if type(S, specfunc(Sum)) then
    return geomsum(sumtoct(op(1, S), v), op(2, S));
  elif type(S, `+`) then
    return map(sumtoct, S, v);
  elif type(S, `*`) then
    L := map(sumtoct, convert(S, list), v);
    return convert(zip(subs, packvars(L, v), L), `*`);
  elif type(S, `^`(anything, posint)) then
    L := [sumtoct(op(1, S), v) $ op(2, S)];
    return convert(zip(subs, packvars(L, v), L), `*`);
  elif type(S, specfunc(CT)) then
    return subs(op(2,S)=v[1], op(1, S));
  elif type(S, specfunc(gfnames)) then
    return eval(gfdict[op(0,S)](v, op(S)));
  else
    return S;
  end if;
end proc;



#### INFINITE SUMS

# Input :
#   - cons, a set of Laurent monomials
#   - G, a directed graph whose vertices are variables and edges are domination
#   relation.
#   - params, set of `small' variables
#
# If possible, returns a directed acyclic graph H extending G such that if the
# variables are ordered according to H, then every monomials in cons is greater
# than 1 (lexicographic ordering).
#
# Raises an error if not possible.
solvecons := proc(cons :: set, G := false, params := {})
  local cons1, rcons, u, v, H;
  uses GraphTheory;

  if G = false then
    H := GraphTheory[Digraph]({seq([1, v], v in indets(cons)), seq(seq([v, u], u in params), v in indets(cons) minus params)});
    return solvecons(cons, H);
  elif not IsAcyclic(G) then error "inconsistent";
  elif nops(cons) = 0 then return G;
  end if;
  
  cons1 := cons[1];
  rcons := cons minus {cons1};
  for v in indets(numer(cons1)) do
    try
      H := CopyGraph(G);
      AddArc(H, {seq([w, v], w in indets(denom(cons1)))});
      return solvecons(rcons, H);
    catch "inconsistent" :
    end;
  end do;

  error "inconsistent";
end proc;

# Returns a linear order compatible with the DAG G.
linearorder := proc(G)
  local racines, H, rest;
  uses GraphTheory;

  if nops(Vertices(G)) = 0 then
    return [];
  end if;

  racines := select(v -> InDegree(G, v) = 0, Vertices(G));
  H := DeleteVertex(G, racines);
  rest := linearorder(H);

  return [op(rest),op(racines)];
end proc;

# Input :
#   - R, a rational function
#   - params, variables assumed to be small
#
# Assume that R is a Laurent formal series w.r.t. variables _W[1],...,_W[r]
# Let T = sum of coefficients of this series
#
# If there exist an order such that T is convergent, then returns T and that order.
# If not, fails
isconvergent := proc(R :: ratpoly, params :: Or(list(name), set(name)))
  local svars, avars, L, ct, co, mord, res, cons, G, ord, den, facts, f; 
  
  svars := indets(R, specindex(_W));
  avars := indets(R) minus svars minus params;
  L := [ op(params), op(avars) ];
  mord :=  plex(op(L));
  cons := {};

  facts := map2(op, 1, factors(denom(normal(R)))[2]);
  for f in facts do
    ct := subs(map(`=`, svars, 0), f);
    co := remove(`=`, [coeffs(collect(normal(f/ct-1), svars, distributed, normal), svars)], 0);
    cons := cons union {seq( Groebner[TrailingTerm](numer(c), mord)[2]/Groebner[TrailingTerm](denom(c), mord)[2], c in co)};
  end do;

  G := solvecons(cons, false, params);
  ord := remove(`=`, linearorder(G), 1);

  mord :=  plex(op(ord));
  res := not `or`(
    seq( Groebner[TestOrder](
      Groebner[TrailingTerm](numer(c), mord)[2],
      Groebner[TrailingTerm](denom(c), mord)[2],
      mord ), c in co ) );

  if res then
    return normal(subs(map(`=`, svars, 1), R)), ord;
  else
    return FAIL;
  end if;
end proc;


# Input :
#   - S, a binomial sum
#   - name, a name
#   - params a set of small parameters
#
# Output :  R, L
#   - R, a rational function
#   - L, a list of variables
#
# The list L gives an order on the variables.
# S is the residue of R with respect to the variables in L that are not
# parameters
sumtores1 := proc(S, name, params)
  local R, vars, x;
  R := sumtoct(S, name);
  vars := select(has, indets(R), name);
  return isconvergent(normal(R/mul(x, x in vars)), params);
end proc;



# Input :
#   - R, a rational function
#   - vars, a list of variables
#   - ord, a positive integer
#
# Returns the first `ord' terms of res_{vars[2..-1]}(R), computed in the Laurent series field
# K((vars[-1]))...((vars[1])).
#
# Very useful to check the consistency of an integral representation.
rser :=
  (R, vars, ord) ->
    map2(foldl, residue, map(normal, series(R, vars[1], ord)), seq(v=0,v in select(member, vars, indets(R))[2..-1])):

# Replaces infinity by maxn in the expression S, replaces Sum by add, Binomial
# by binomial, etc and evaluates.
computesum := (S, maxn) -> eval(subs([Sum=add, infinity=maxn, op(op(op(gfdictcomp)))], S));




###### GEOMETRIC REDUCTION OF PERIODS

inorout_base := proc(S, i)
	local minx, maxx, miny, left, right, dom, middle, ret;
	minx := min(map2(op, 1, S)); 
	maxx := max(map2(op, 1, S));

	if minx=maxx then return {}; end if;
 	if nops(S[1]) <= 1 or i <= 1 then return {-1} end if;

 	miny := min(map2(op, 2, S));
 	middle := map2(subsop, 2=NULL, select(m -> m[2] = miny, S));

 	left := min(map2(op, 1, middle));
  	right := max(map2(op, 1, middle));
	ret := inorout_base(middle, i-1);

	if left > minx then ret := ret union {1}; end if;
	if right < maxx then ret := ret union {-1}; end if;

	return ret;
end proc;


inorout := proc(P, T, ord)
  local S, vars;
  vars := remove(has, ord, T);
  coeffs(collect(P, [T, op(vars)], distributed), [T, op(vars)], 'mon');
  S := map2(map2, degree, {mon}, [T, op(vars)]);

  return inorout_base(S, ListTools[Search](T, ord));
end proc;


# Hermite reduction
# Input :
#   - R, a rational function
#   - v, a name
#   - cert, a boolean (default: false)
#
# Output :
#   A rational function S such that R - S = T' for some rational function T'
#   and such that S has poles of order at most 1 w.r.t. v (including at infinity)
#
#   If cert=true, then it also returns T.
hermitered := proc(R, v :: name, cert := false)
  local a, d, g, dm, ds, dm2, dms, b, c, k;
  a := numer(R);
  d := denom(R);
  
  g := 0;
  dm := gcd(d, diff(d, v));
  dm := normal(dm);
  ds := normal(d/dm);

  while degree(dm, v) > 0 do
    dm2 := gcd(dm, diff(dm, v));
    dms := normal(dm/dm2);
    gcdex(-normal(ds*diff(dm,v)/dm), dms, a, v, 'b', 'c');
    a := normal(c - diff(b,v)*ds/dms);
    g := g + b/dm;
    dm := dm2;
  end do;
  
  if cert then
    return normal(a/ds/dm), g;
  else
    return normal(a/ds/dm);
  end if;
end proc;

# Input:
#   - R, rational function
#   - v, symbol
#   - ord, list of symbols containing v
#
# Returns FAIL or a rational function without v which is the sum of the
# residues of R at `small' poles.
#
# [Implementation using hermitered and Rothstein-Tragger resutant
#   for the residue computation]
ratres := proc(R, v, ord)
  local Rn, F, tot, f, ioo, Q, res, n;

  if type(R, `+`) then
    return map(ratres, R, v, ord);
  end if;

  Rn := hermitered(normal(R), v);
  F := map2(op, 1, select(has, factors(denom(Rn))[2], v));
  
  tot := 0;
  for f in F do
    if type(f, `+`) then
      ioo := inorout(f, v, ord);

      if ioo = {1} then
        Q := normal(Rn*f/diff(f, v));
        res := collect(resultant(numer(Q)-Z_*denom(Q), f, v), Z_);
        n := degree(res, Z_);
        tot := tot - normal(coeff(res, Z_, n-1)/coeff(res, Z_, n));
      elif nops(ioo) > 1 then
        return FAIL;
      end if;
    else  # f = v
      tot := tot + normal(subs(v=0,normal(v*Rn)));
    end if;
  end do;

  return normal(tot);
end proc;

# Idem
#
# [Implementation using maple's residue(R, v=infinity)]
ratres2 := proc(R, v, ord)
  local terms, F, tot, f, ioo, Q, res, n;

  terms := convert(R, parfrac, v);
  if type(terms, `+`) then
    terms := [op(terms)];
  else
    terms := [terms];
  end if;

  tot := 0;
  for F in terms do
    f := denom(F);
    if eval(f, v=0) = 0 then
      tot := tot + residue(F, v=0);
    elif degree(f, v) > 0 then
      ioo := inorout(f, v, ord);
      if ioo = {1} then
        tot := tot - residue(F, v=infinity);
      elif nops(ioo) > 1 then
        return FAIL;
      end if;
    end if;
  end do;

  return normal(tot);
end proc;


geomred := proc(R :: ratpoly, ord :: list, params :: set)
  local S, red, v;
  S := R;
  for v in indets(R) minus params do
    red := ratres(S, v, ord);
    if red <> FAIL then S := red; end if;
  end do;
  return collect(factor(S), params, factor, distributed);
end proc;

geomredall := proc(R :: ratpoly, ord :: list, params :: set)
  local tbl, cur, ivars, st, red, v, nst, all, deg;

  tbl := table([ {} = normal(R) ]);
  cur := { {} };
  ivars := indets(R) minus params;

  while nops(cur) > 0 do
    st := cur[1];
    red := tbl[st];
    if red <> FAIL then
      for v in ivars minus st do
        nst := st union {v};
        if type(tbl[nst], indexed) then
          tbl[nst] := factor(ratres(red, v, ord));
        end if;
        cur := cur union {nst};
      end do;
    end if;
    cur := cur minus {st};
  end do;
  
  all := remove(`=`, map(op, [entries(tbl)]), FAIL);
  cur := min(map(nops@indets, all));
  all := select(f -> nops(indets(f))=cur, all);
  deg := f -> degree(denom(f), ivars) + min(0, degree(numer(f),ivars)-degree(denom(f),ivars) +nops(ivars)+1);
  cur := min(map(deg, all));
  all := select(f -> deg(f)=cur, all);

  return map(f -> collect(factor(f), params, factor, distributed), all);
end proc;

# Input :
#   - S, a binomial sum
#   - name, a name
#   - params a set of small parameters
#
# Output :  R, L
#   - R, a rational function
#   - L, a list of variables
#
# The list L gives an order on the variables.
# S is the residue of R with respect to the variables in L that are not
# parameters
sumtores := proc(S, name, params)
  local R, ord;

  R, ord := sumtores1(S, name, params);
  R := geomred(R, ord, params);
  ord := select(has, ord, indets(R));

  return op(subs(packvars([ord], name)[1], [R,ord]));
end proc;


end module:

