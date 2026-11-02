# Literal Appendix A.1 FRLP: backward feasibility and safe pruning

The target LP is unchanged: each rate leaf uses its lower endpoint in A.2, its upper endpoint in A.3, and ordinary triangular interpolation at the exact A.2 state. Numerical formulas below are only lower bounds for those fixed rows, not replacement LP rows.

Let P be the ordinary PL interpolant of the already fixed next-layer table, with all diagonal values zero and all entries nonpositive. Enforce numerical coordinate monotonicity by lowering entries only. The actual source uses these two exact integer passes:

```
for i=1,...,m:
    for j=i-1,...,0:
        v[i,j] = min(v[i,j], v[i,j+1])
for i=1,...,m:
    for j=0,...,i-1:
        v[i,j] = min(v[i,j], v[i-1,j])
```

Each pass preserves every current LP row because its RHS uses only the already fixed next layer and each changed left-side value decreases. The first pass makes rows nondecreasing. The second pass makes columns nonincreasing and preserves row monotonicity: each new row is the pointwise minimum of two nondecreasing rows (with the diagonal value zero). Thus the output is nondecreasing in x and nonincreasing in y. On each of the two grid triangles, the corresponding partial derivatives of P are differences of neighboring entries divided by positive knot widths, so P has the same coordinate monotonicity.

## A lower bound for all literal leaf continuations in a tree node

Write

- E(r)=exp(-Delta r), I(r)=(1-E(r))/r, I(0)=Delta;
- Phi(r)=1-y+y[t+I(r)+E(r)(1-t-Delta)];
- Y(r)=y E(r)/Phi(r);
- R(r)=x/y+b Delta (r+exp(-r)-1)/c;
- X(r)=Y(r) min(1,R(r)).

Phi and Y are nonincreasing in r. For Y, Phi/E = (1-y)exp(Delta r)+y[t exp(Delta r)+I(r)exp(Delta r)+1-t-Delta], whose terms are nondecreasing; I(r)exp(Delta r)=integral_0^Delta exp(ru)du. R is nondecreasing.

For r in [a,b], if R(a)>=1 then all continuations are zero. Otherwise set Y_b=Y(b), or zero when b=infinity, and X_low=Y_b R(a). Then X(r)>=X_low even after clipping, and Y(r)<=Y(a). Hence, using coordinate monotonicity and P<=0,

Phi(r) P(Y(r),X(r)) >= Phi(a) P(Y(a),X_low).

This lower bound applies to all *literal leaf* rows inside the tree node; it does not select a different ancestor affine row.

## Running-cost lower bound across a tree node

For the leaf upper endpoint r, C(r)=min(p Delta,p I(r)), alpha(r)=y b E(r) Delta^2+(-p)^+ I(r) Delta, beta(r)=Delta[b(1-y+ty)+y b I(r)] are all nonincreasing in r.

The piecewise free bound H0 is beta-beta^2/(4 alpha) for beta<=alpha, and alpha[3/4+(1/2)log(beta/alpha)] for beta>alpha. It is increasing in each positive argument. The polynomial version Hpoly replaces the log part by alpha f(w), f(w)=3/4+w+w^3/3+w^5/5 and w=(beta-alpha)/(beta+alpha). It too is coordinate-increasing: in the log branch

partial_alpha Hpoly = 1/4+w+w^3/3+w^5/5+w^6/2 > 0,
partial_beta Hpoly = (1-w)^2(1+w^2+w^4)/2 >= 0.

Thus C(b)+H0(alpha(b),beta(b)), or the weaker polynomial version, is below every descendant leaf's full A.3 cost. The shifted bound itself need not be monotone. It can nevertheless be enclosed uniformly on a node without changing any leaf row. Put z(r)=beta(r) exp(-r), also nonincreasing. If z(b)>alpha(a), every descendant leaf satisfies the original shifted guard, and its shifted bound is at least

C(b)+z(b)+alpha(b) a/2-alpha(a)^2/[24(z(b)-alpha(a)/2)].

This follows by bounding each positive term below and the negative fraction above; its denominator is positive. Use this only when b<=64 and the conservative guard has eps slack. With independent variables A=alpha(a), B=alpha(b), Z=z(b), the derivatives of F=Z+B a/2-A^2/[24(Z-A/2)] obey partial_Z F<=7/6, |partial_A F|<=1/4, partial_B F=a/2<=32 under Z>A. Thus input and operation errors are below400 eta, well covered by eps; no monotonicity of the shifted formula is assumed. The solver may additionally use alpha(a)+16 eta, max(0,alpha(b)-16 eta), max(0,z(b)-16 eta), and max(0,r_a-16 eta) as conservative inputs. At infinity use min(p Delta,0).

## Fast binary64 error allowance

Assume validated frozen input errors below eta=2^-48, relative error of Delta below eta, 1/4<b<2, .85<c<.9, finite rates at most8192, and strict round-to-nearest binary64 with no contraction/reassociation. Let eps=2^-32=65536 eta, U=2^40.

Compute Rhat=x/y+bhat Deltahat vhat(a)/chat. Since x/y<=1, b Delta v<=25.6, c>.85, direct propagation gives |Rhat-R|<512 eta (a loose bound; the elementary estimate is below100 eta). Therefore Rhat>=1+eps safely implies R>=1. Otherwise R<1+eps+512 eta and b Delta(a-1)/c<=R implies a Delta<4.

For hi, use Y_b=0 unless the test rhat_hi Deltahat<=4-eps holds. This test implies the exact hi Delta<4; otherwise zero is always a safe lower choice. In retained divisions, E>exp(-4)>1/64 and Phi>=E. Existing frozen-coefficient propagation gives |Phihat-Phi|<32 eta and |Yhat-Y|<4096 eta for either endpoint. The product Yhat_b Rhat has error below8192 eta. The original interpolation routine's shifts y+eps and x-eps therefore place the integer query above the exact Y(a) and below the exact X_low, with ample slack. Its exact integer barycentric floor gives a lower bound for P at that query. The original upper multiplier ceil(U(Phihat+eps))+1 bounds U Phi(a) from above; multiplication by the nonpositive interpolation preserves a lower bound.

The original cost arithmetic allowance eps already dominates1024 eta. Adding the full free logarithmic bound is safe with the existing alpha,beta>eps guard: exact beta<.004 and alpha>eps-8 eta imply beta/alpha<2e7. Hence partial_alpha H0<9 and partial_beta H0<=1; coefficient perturbations contribute below100 eta. A range-reduced positive atanh expansion with12 terms is a deliberate lower approximation to log. With ratio>1, frexp reduces its mantissa to[1,2), and the omitted tail is positive. Use ln2 lower constant0x1.62e42fefa39eep-1, independently verified below ln2 using40 rational atanh terms. The exponent is at most24. Basic operation errors in this finite log formula are well below1024 eta; subtracting eps from the computed log and retaining the original eps cost allowances covers all upward errors. No platform log accuracy assumption is needed.

## Backward feasibility

At each state the source starts with the minimum of zero and conservative leaf scores. A leaf score is its computed conservative RHS minus an additional exact integer slack

    delta = 2^-28 = 4096 / 2^40.

The selected source also subtracts a nonnegative integer allowance from every proposed or visited leaf score. It sets w=y-x, computes d=clamp(w-E(r)*max(0,w-y*b*Delta*v(r)/c),0,w), and subtracts ceil(2^24*d)/U. The clamp gives 0<=d<=w<=1, so this additional subtraction is nonnegative and at most epsilon+1/U with epsilon=2^-16. Therefore every resulting candidate is still no larger than its original leaf-row RHS. This expression is only a policy for proposing conservative values: it does not enter any original row or node pruning bound, and no accuracy or telescoping claim about d is used in the proof.

The two initial leaves are the previous time-layer hint and the unbounded last leaf. If the root cannot yet be pruned, eighteen ternary-search proposal iterations and three remaining interval locations supply additional candidate leaf scores. These are only proposals: the final `search(0,Q,...)` still covers the entire original leaf set. No unimodality assumption or ternary-search correctness assumption enters the proof.

The recursive search takes minima with conservative evaluated leaf RHS values minus delta and the nonnegative allowance, and prunes a node only when its genuine lower bound for every literal descendant row is at least the current best. Visited leaves and pruned nodes cover all leaves including the unbounded one. Therefore the final raw value is no larger than every literal LP row RHS. Subsequent monotonicity lowering preserves these inequalities. Starting from terminal zeros proves full LP feasibility by backward induction. This is not a claim that monotonicity raising preserves LP feasibility.


## Exact finite coverage, integer range, and execution conditions

Each recursive split partitions [lo,hi) into two disjoint child index intervals; the final unbounded leaf is included in [0,Q). When a node is pruned the source adds hi-lo to its coverage counter; when a leaf is evaluated it adds one. At every state `coverage==Q` is asserted. The recursive structure shows this is a disjoint cover, and the proof above establishes that each covered row has RHS at least the final best. Earlier best candidates can only decrease afterward, so pruning remains valid. There is no maximum over ancestor affine rows used as a leaf score.

Let M_k be the largest absolute table value at layer k after division by U=2^40. The baseline cost, the continuation multiplier, downward integer floors, delta, and the allowance bound epsilon+1/U give

M_k <= (1+2 eps+2/U) M_(k+1) + 2/640 + 2 eps + 5/U + delta + epsilon,    M_640=0.

The monotonicity-lowering passes do not change the global minimum below the pre-pass global minimum. Therefore this recurrence bounds them too. With delta=2^-28 and epsilon=2^-16, its geometric-series bound is below2.009768912, hence below2.01 and certainly below4. In particular every stored table value has magnitude below4U. Grid widths and offsets are at most2^30; barycentric products fit signed64, and weighted table numerators have magnitude at most4U*2^60=2^102, below signed128 capacity. Continuation multiplier-table products are below2^84. The exact negative-floor implementation and array sizes therefore remain safe.

The source rejects fast-math and NDEBUG. It sets and checks FE_TONEAREST on every OpenMP outer-loop work item, hence on each worker before the state computations. Compile with contraction disabled (`-ffp-contract=off`) and ordinary strict binary64 operations. The frozen inputs and their rational coefficient checks are those of the existing validated header and original rational threshold/rate/grid files.

The result numerator is HU+v[0,m,0], where HU=floor(U theta(0)); therefore gamma=(HU+v[0,m,0])/U <= theta(0)+v[0,m,0]. Together with all row inequalities, terminal zeros, diagonal zeros, and nonpositivity, this is a feasible solution of the original fixed-coefficient FRLP.

Audited source: `iid_literal_solver/literal_potential16.cpp`, SHA-256 `6cafb1fb9e1be830c8b2c8cd60304bb1e2c3496e3da5a3e7e2d47cb4b9a84219`.

The exact rate file contains the endpoint64, represented exactly in binary64, and no other endpoint within10^-6 of64. Thus the `rr<=64` scalar guard cannot admit an exact endpoint above64 through input rounding. The last finite endpoint is exactly8192. These fixed-data facts are checked in `audit/ARITHMETIC_CHECKS.json`.
