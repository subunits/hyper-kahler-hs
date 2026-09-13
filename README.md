# AVS Hyper-Kahler Extension

Single-file Haskell playground implementing the quaternionic (H^k) layer
of the Augmented Vector Space stack.

Paste `HyperKahler.hs` into https://play.haskell.org and hit Run.
Requires only `base`. No cabal, no dependencies.

---

## Stack position

```
avsp            R^n, L2, lag augmentation, kNN, OLS, entropy, anomaly
  kahler-ts     C^(n/2), Hermitian H=g+iw, Chern proxy, VR homology
    kahler_isco Kerr ISCO geodesic coords as Kahler signals
      bh-phase-space  MP horizons, black rings, GL instability, ZAMO
        HyperKahler.hs  [THIS FILE]
```

Each layer strictly extends the one below.  The Riemannian metric g from
`avsp` is preserved exactly as Re H throughout — no information is lost
in the complexification or quaternification steps.

---

## Key idea

Each ISCO time step (r, t, theta, phi) fills one quaternion natively:

```
Q r t theta phi
```

No lag augmentation needed.  The four spacetime coordinates are already
the right dimension for H^1.  The three symplectic forms that fall out of
the quaternionic Hermitian inner product correspond directly to the three
Killing-Yano tensors of Myers-Perry geometry (Frolov-Kubiznak 2007):

```
H(u,v) = g(u,v) + i*omegaI(u,v) + j*omegaJ(u,v) + k*omegaK(u,v)

omegaI  ~  Omega_I  (radial Killing-Yano tensor)
omegaJ  ~  Omega_J  (polar  Killing-Yano tensor)
omegaK  ~  Omega_K  (azimuthal Killing-Yano tensor)
```

---

## What it computes

### Quaternion algebra verification
I^2 = J^2 = K^2 = IJK = -1 checked on four test quaternions.
Non-commutativity: i*j = +k, j*i = -k confirmed.

### ISCO native embedding
Each of the 60 ISCO time steps embedded as Q(r, t, theta, phi).
Quaternionic norm |q| reported alongside the four coordinates.
Re H(v,v) = ||v||^2 confirms the Riemannian metric is preserved.

### Hyper-Kahler conditions
Checked on consecutive pairs of embedded points:
- wI skew-symmetric:   wI(u,v) = -wI(v,u)
- wI compatible with I: wI(u^,v^) = g(I u^, v^)  (on normalised vectors)
- same for wJ and wK

Compatibility is checked on unit-normalised vectors to separate two
distinct failure modes:
- FAIL after normalisation  = genuinely non-Lagrangian embedding
- FAIL only before normalise = coordinate scale artefact

Raw residuals resI, resJ, resK are reported for every pair.

### Curvature d(omega) over triples
Discrete exterior derivative of each symplectic form computed over
consecutive triples of points.  Non-zero = curved HK manifold.
ISCO shows monotonically growing dωJ and dωK — curvature accumulates
as the geodesic winds.

### Chern proxies
Discrete integral sum_{i<j} w(p_i, p_j) for each of the three
symplectic forms, computed per signal and displayed as bar charts.

Signals compared:
- ISCO geodesic (60 points, native H^1)
- MP d=5 horizon data (4 points, H^1)
- Black ring data (5 points, H^1)
- MP d=6 single-spin baseline (4 points, H^2)
- MP d=6 multi-spin regimes A/B/C (4 points each, H^2)

### VR persistent homology in H^k
Vietoris-Rips filtration in Hermitian L2 distance.
Reports H0 (connected components) and H1 (loops) at 20 epsilon steps.
All 60 ISCO points used.  H1=0 throughout because the data covers
less than one full orbital period.

### MP d=6 multi-spin (three regimes)
Two quaternions per point in H^2:
  [Q(r_H, a1, a2, A_H),  Q(mu, nu1, nu2, Omega)]

Three spin regimes tested:
- A: a1 >> a2  (near single-spin)
- B: a1 = a2   (equal spins, SU(2) symmetry)
- C: a2 = 0    (single-spin limit, recovers mp6 baseline exactly)

### Spin symmetry fingerprint
Classifies each signal by its Chern proxy ratios:
- single-spin (I only):  cJ = cK = 0
- near SU(2) [cK=0]:     cI ~ cJ, cK = 0  (equal-spin case)
- near SU(2):            cI ~ cJ ~ cK
- I-dominant:            |cI| > |cJ|, |cI| > |cK|
- mixed:                 none of the above

---

## Results summary

| Signal             | omegaI    | omegaJ   | omegaK    | Reading              |
|--------------------|-----------|----------|-----------|----------------------|
| ISCO               | -1161     | -52      | -7227     | azimuthal dominant   |
| MP d=5             | 1.5       | -0.3     | 53.8      | weak azimuthal       |
| Black ring         | 903       | 515      | -899      | I/K competing        |
| MP d=6 1-spin      | 132       | 0.0      | 0.0       | I only (Lagrangian)  |
| MP d=6 A (a1>>a2)  | 7.2       | -2.3     | 2.6       | I-dominant           |
| MP d=6 B (a1=a2)   | 4.96      | 4.96     | 0.0       | near SU(2) [cK=0]   |
| MP d=6 C (a2=0)    | 132       | 0.0      | 0.0       | single-spin, exact   |

Key findings:
- ISCO is non-Lagrangian in H^k: resI ~ -resK (I/K antisymmetry) throughout
- MP d=6 single-spin is Lagrangian in the I complex structure; cJ=cK=0 exactly
- MP d=6 equal-spin shows dωI = dωJ, dωK = 0 at every triple — SU(2) in curvature
- Regime C recovers the single-spin baseline to 4 decimal places — strongest
  internal consistency check in the stack
- Black ring: large competing I and K Chern proxies reflect the two rotation
  directions (S1 ring angle and S2 sphere angle) pulling against each other

---

## HK conditions table

| Signal          | wI-skew | wI-compat | wJ-compat | wK-compat | geometry           |
|-----------------|---------|-----------|-----------|-----------|--------------------|
| ISCO            | OK      | FAIL      | FAIL      | FAIL      | non-Lagrangian     |
| MP d=5          | OK      | FAIL      | mixed     | FAIL      | non-Lagrangian     |
| Black ring      | OK      | FAIL      | FAIL      | FAIL      | non-Lagrangian     |
| MP d=6 1-spin   | OK      | OK        | OK        | OK        | Lagrangian (I)     |
| MP d=6 A        | OK      | OK        | OK        | FAIL      | near-Lagrangian    |
| MP d=6 B        | OK      | OK        | OK        | OK        | Lagrangian (I+J)   |
| MP d=6 C        | OK      | OK        | OK        | OK        | Lagrangian (I)     |

---

## Known limitations

- VR H1 = 0 throughout: ISCO data covers < 1 full orbit.  Need ~120 points
  with full phi closure to detect the H1 loop generator.
- MP horizon data is sparse (4 points).  Chern proxy magnitudes are
  meaningful for comparison but not absolute.
- Multi-spin d=6 horizon radii are approximated.  For exact values, feed
  output from bh-phase-space bisection directly.
- Spin symmetry fingerprint uses a divide-by-zero guard (999) when cK=0.
  Replace with direct threshold classification for production use.
- No GL instability check in the HK layer yet.  The membrane limit
  (ultra-spinning MP) should show a sharp change in Chern proxy ratios.

---

## Parameter conventions

| Symbol   | Meaning                                           |
|----------|---------------------------------------------------|
| Q        | Quaternion q0 + q1*i + q2*j + q3*k               |
| H^k      | Quaternionic vector space, k quaternions per point |
| H(u,v)   | Quaternionic Hermitian inner product               |
| g        | Re H  (Riemannian metric, = L2 from avsp)         |
| omegaI   | q1 component of H  (Killing-Yano Omega_I)         |
| omegaJ   | q2 component of H  (Killing-Yano Omega_J)         |
| omegaK   | q3 component of H  (Killing-Yano Omega_K)         |
| nu       | r_H / a  (parametric spin variable)               |
| cI/cJ/cK | Chern proxies: sum_{i<j} omegaX(p_i, p_j)        |
| resI/J/K | Compatibility residuals on normalised vectors     |

---

## What comes next

1. **Full ISCO orbit** — ~120 points with phi returning to start.
   Prediction: H1 generator appears in VR at epsilon ~ orbital diameter.

2. **d=5 ring vs MP triple-degeneracy VR comparison** — feed the a_H(j)
   curve points for thin ring, plump ring, and MP branches as separate
   QVec signals.  Prediction: topologically distinct H0/H1 structure
   across the j_cusp boundary.

3. **Spin symmetry fingerprint v2** — replace ratio-based classification
   with direct threshold on |cK| and |cI - cJ| to eliminate the 999
   divide-by-zero artefact and give clean three-way discrimination.

4. **GL instability in H^k** — embed the (r_H, a, k, Im_Omega) tuple
   as a quaternion and track how the Chern proxies change as k crosses
   k_GL.  Prediction: sharp omegaK transition at onset.

---

## Source references

1. R. C. Myers and M. J. Perry, Ann. Phys. 172, 304 (1986)
2. R. Emparan and H. S. Reall, Phys. Rev. D 65, 084025 (2006)
3. V. P. Frolov and D. Kubiznak, Phys. Rev. Lett. 98, 011101 (2007)
4. D. Pereñiguez Rodriguez, arXiv:1808.04009 (2018)
5. avsp -- github.com/subunits/avsp
6. kahler-ts -- github.com/subunits/kahler-ts
7. bh-phase-space -- github.com/subunits/bh-phase-space
