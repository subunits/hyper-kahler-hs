# hyper-kahler-hs

Augmented Vector Space — Hyper-Kähler Extension. The top layer of a five-repo
stack that lifts signal analysis from R^n through C^(n/2) to H^(n/4), where
the three symplectic forms of the quaternionic Hermitian metric correspond to
the Killing-Yano tensors of Myers-Perry black hole geometry.

Runs on [play.haskell.org](https://play.haskell.org) with no dependencies
(base only). Compile locally with `ghc -O HyperKahler.hs -o hk`.

---

## Stack

```
avsp  →  kahler-ts  →  kahler_isco  →  bh-phase-space  →  hyper-kahler-hs
R^n      C^(n/2)       C^(n/2)          physics              H^(n/4)
```

Each layer extends the one below:

| Repo | Space | Key addition |
|---|---|---|
| `avsp` | R^n | Lag augmentation, kNN, entropy |
| `kahler-ts` | C^(n/2) | Hermitian H=g+iω, Chern proxy, VR homology |
| `kahler_isco` | C^(n/2) | ISCO geodesic coords as Kähler signal |
| `bh-phase-space` | — | Myers-Perry, black rings, GL instability, ZAMO |
| `hyper-kahler-hs` | H^(n/4) | Three symplectic forms, Killing-Yano correspondence |

---

## Key insight

Each ISCO time step `(r, t, θ, φ)` fills one quaternion natively — no lag
augmentation needed. The Hermitian inner product

```
H(u,v) = Σ conj(u_a) · v_a  ∈  H
```

splits into four real parts:

```
H = g + i·ωI + j·ωJ + k·ωK
  g   = Re H   Riemannian metric  (= original L2 from avsp)
  ωI  = q1 H   Symplectic form I  (~ Killing-Yano Ω_I, radial)
  ωJ  = q2 H   Symplectic form J  (~ Killing-Yano Ω_J, polar)
  ωK  = q3 H   Symplectic form K  (~ Killing-Yano Ω_K, azimuthal)
```

The three complex structures I, J, K are left-multiplication by i, j, k and
satisfy I²=J²=K²=IJK=-1 (verified at runtime on arbitrary quaternions).

Reference: Frolov & Kubiznak (2007), hidden symmetries of Myers-Perry geometry.

---

## Physics

### Horizon solvers

All horizon radii are computed from first principles; nothing is hardcoded.
Change `bigG` or `massM` at the top of the file and every horizon, area, and
angular velocity updates automatically.

**d=5 single-spin** — exact closed form:
```
r_H = sqrt(μ - a²)
```

**d=6 single-spin** — Cardano solution of the depressed cubic `r³ + a²r - μ = 0`:
```
r_H = cbrt(-q/2 + √disc) + cbrt(-q/2 - √disc)
      where q = -μ,  disc = (q/2)² + (a²/3)³
```

**d=6 two-spin** — closed-form quadratic in u = r²:
```
u² + (a₁² + a₂² - μ)u + a₁²a₂² = 0
r_H = sqrt(larger non-negative root)
Extremal bound: a₁ + a₂ ≤ sqrt(μ) ≈ 0.691  (M=G=1)
```

Mass parameter:
```
μ = 16π G M / [(d-2) Ω_{d-2}]
```

where Ω_n = 2π^{(n+1)/2} / Γ((n+1)/2) uses a clean two-base-case recursion
(same fix as bh-phase-space).

### Signals embedded as quaternion points

| Signal | Embedding | Notes |
|---|---|---|
| ISCO geodesic | Q(r, t, θ, φ) per step | Native — no augmentation |
| MP d=5 horizon | Q(r_H, a, A_H, Ω_H) | Computed from μ₅ |
| Black ring | Q(ν, λ, A_H, Ω_H) | λ = equilibrium value |
| MP d=6 single-spin | [Q(r_H, a, 0, 0), Q(μ, ν, 0, 0)] | H² embedding |
| MP d=6 two-spin A | [Q(r_H, a₁, a₂, 0), Q(μ, ν₁, ν₂, 0)] | a₁ >> a₂ |
| MP d=6 two-spin B | same | a₁ = a₂ (equal spins) |
| MP d=6 two-spin C | same | a₂ = 0 (single-spin limit) |

### Hyper-Kähler condition checks

Compatibility is checked on unit-normalised vectors to separate genuine
non-Lagrangian embeddings from coordinate scale artefacts:

- **FAIL after normalisation** = genuinely non-Lagrangian
- **FAIL only before** = coordinate scale mismatch (physical units)

ISCO FAILs with resI ≈ −resK (I/K antisymmetry) — the orbit is a curved
embedding in H^k, not a Lagrangian submanifold.

MP d=6 single-spin and equal-spin (regime B) pass all conditions with
residuals exactly zero.

### Multi-spin regimes (d=6, M=G=1)

All spin values are within the extremal bound a₁+a₂ ≤ sqrt(μ₆) ≈ 0.691.

| Regime | Spins | Chern result | Interpretation |
|---|---|---|---|
| A (a₁ >> a₂) | (0.60, 0.05) … | cK dominant | Azimuthal KY from large a₁ |
| B (a₁ = a₂) | (0.30, 0.30) … | cI = cJ, cK = 0 | SU(2) symmetry |
| C (a₂ = 0) | (0.0, 0.3, 0.5, 0.6) | cJ = cK = 0 exactly | Single-spin limit |

Regime B equality cI=cJ with cK=0 is the SU(2) fingerprint: the two equal
rotation planes are exchanged by the symmetry, and the K form that
distinguishes them vanishes.

---

## VR persistent homology

Uses an incremental Vietoris-Rips filtration in the Hermitian L² distance on
H^k:

1. Distance matrix built once as `Array (Int,Int) Double` — O(1) lookup
2. All edges sorted once — O(n² log n)
3. Epsilon steps walk the sorted list incrementally; state (union-find, edge
   array, triangle count) is carried forward, never rebuilt
4. Triangle count: for new edge (i,j), count k with both (i,k) and (j,k)
   already present — O(n) per new edge

Cost: O(n² log n + n³) spread across the filtration, vs the original
O(steps × n³) which repeated full triangle enumeration at every epsilon.

**ISCO result**: H0 drops from 59→1 by ε≈1.98 (all 60 points merge into one
component). H1 spikes at ε≈0.84 (36 loops) then collapses to 0. The orbit
covers less than one full revolution in the sampled data — denser sampling is
needed to recover the H1 generator from the closed periodic orbit.

---

## Bugs fixed in this version

1. **`gammaH` in sphere area** — replaced fragile factorial recursion with
   clean two-base-case recursion from Γ(1/2)=√π and Γ(1)=1 (inherited fix
   from bh-phase-space).

2. **`vrFiltration` O(steps × n³)** — replaced with incremental walk using
   `Array`-backed union-find and edge membership for O(1) lookup.

3. **`zamoOmega` naming** — `zamoOmega5D` renamed `zamoOmega5DEqual` with
   comment clarifying it assumes a=b; discrepancy near horizon is expected, not
   a bug (fix in bh-phase-space).

4. **Hardcoded horizon radii** — `mpPoints`, `mp6Points`, `mpA`, `mpB`, `mpC`
   previously had hardcoded `r_H` values copied from Playground.hs output.
   Now computed at runtime from `bigG`/`massM` via exact closed-form solvers.
   The old `mpA`/`mpB` spin values also exceeded the extremal bound (no horizon
   exists for a₁+a₂ > sqrt(μ)) and have been replaced with physically valid
   values.

---

## Next steps

- **ISCO full orbit** — denser sampling (≥60 steps per revolution) to recover
  the H1 generator; the periodic orbit is a loop in H^k that should produce a
  persistent 1-cycle.
- **d=5 ring vs MP triple-degeneracy** — VR comparison near j=j_cusp where
  the thin ring, plump ring, and MP branches all meet; the topology of the
  phase space should be visible in H0/H1.
- **Regime A physics** — the K-dominance in the a₁>>a₂ regime warrants deeper
  analysis; the azimuthal KY tensor should be suppressed, not enhanced, when
  one spin vanishes. May indicate the quaternion slot assignment (which
  coordinate maps to which q component) needs revisiting against the
  Frolov-Kubiznak tensor conventions.
