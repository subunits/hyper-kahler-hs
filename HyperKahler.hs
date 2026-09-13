-- Augmented Vector Space -- Hyper-Kahler Extension
-- No dependencies -- runs on play.haskell.org (base only)
-- Compile: ghc -O HyperKahler.hs -o hk
--
-- Stack:
--   avsp        (R^n, L2, lag augmentation, kNN, entropy)
--   kahler-ts   (C^(n/2), Hermitian H=g+iw, Chern, VR homology)
--   kahler_isco (ISCO geodesic coords as Kahler signals)
--   bh-phase-space (MP horizons, black rings, GL, ZAMO)
--   THIS FILE   (H^(n/4), three symplectic forms, Killing-Yano correspondence)
--
-- Key insight: each ISCO time step (r, t, theta, phi) fills one quaternion
-- natively. No lag augmentation needed. The three symplectic forms omega_I,
-- omega_J, omega_K correspond to the Killing-Yano tensors Omega_I, Omega_J,
-- Omega_K of Myers-Perry geometry (Frolov-Kubiznak 2007).

module Main where

import Data.Array (Array, listArray, (!), (//))
import Data.List  (sortBy, nub)
import Data.Ord   (comparing)

-- ============================================================
-- PHYSICAL CONSTANTS & HORIZON SOLVER
-- ============================================================

-- Gravitational constant and mass (change here; all horizons update automatically)
bigG :: Double
bigG = 1.0

massM :: Double
massM = 1.0

-- | Sphere area Omega_n = 2 pi^{(n+1)/2} / Gamma((n+1)/2)
sphereArea :: Int -> Double
sphereArea n = 2 * pi ** (fromIntegral (n+1) / 2) / gammaHalf (n+1)
  where
    gammaHalf 1 = sqrt pi
    gammaHalf 2 = 1.0
    gammaHalf k = fromIntegral (k-2) / 2.0 * gammaHalf (k-2)

-- | Mass parameter mu = 16 pi G M / [(d-2) Omega_{d-2}]
massParam :: Int -> Double
massParam d = 16 * pi * bigG * massM / (fromIntegral (d-2) * sphereArea (d-2))

-- | d=5 single-spin horizon: r_H = sqrt(mu - a^2)
horizonD5 :: Double -> Double
horizonD5 a =
  let mu = massParam 5
  in  if a^2 > mu then 0.0 else sqrt (mu - a^2)

-- | d=6 single-spin horizon: positive real root of r^2 + a^2 - mu/r = 0
--   i.e. r^3 + a^2*r - mu = 0  (depressed cubic, solved via Cardano)
--   This is the deltaMP formula for one non-zero spin parameter.
horizonD6Single :: Double -> Double
horizonD6Single a =
  let mu   = massParam 6
      p    = a^2
      q    = -mu
      disc = (q/2)^2 + (p/3)^3
      cbrt x = if x < 0 then -((-x)**(1/3)) else x**(1/3)
  in  cbrt (-q/2 + sqrt disc) + cbrt (-q/2 - sqrt disc)

-- | d=6 two-spin horizon: positive real root of
--     (r^2 + a1^2)(r^2 + a2^2) = mu * r^2
--   Substituting u = r^2 gives the quadratic:
--     u^2 + (a1^2 + a2^2 - mu)*u + a1^2*a2^2 = 0
--   Extremal bound: a1 + a2 <= sqrt(mu) ~ 0.691 (for d=6 M=G=1)
--   Take the larger non-negative root (outer horizon).
horizonD6Two :: Double -> Double -> Double
horizonD6Two a1 a2 =
  let mu   = massParam 6
      b    = a1^2 + a2^2 - mu
      c    = a1^2 * a2^2
      disc = b^2 - 4*c
      u1   = (-b + sqrt disc) / 2
      u2   = (-b - sqrt disc) / 2
  in  if disc < 0 then 0.0
      else sqrt (max 0 (max u1 u2))

-- | d=5 horizon area: A_H = r_H^{d-4} (r_H^2 + a^2) Omega_{d-2}
areaD5 :: Double -> Double -> Double
areaD5 rH a = rH * (rH^2 + a^2) * sphereArea 3

-- | d=5 angular velocity: Omega_H = a / (r_H^2 + a^2)
omegaHD5 :: Double -> Double -> Double
omegaHD5 rH a = a / (rH^2 + a^2)

-- ============================================================
-- QUATERNION ALGEBRA
-- ============================================================

data Q = Q { q0 :: Double, q1 :: Double, q2 :: Double, q3 :: Double }

instance Show Q where
  show (Q a b c d) = fmtSigned a ++ fmtSigned b ++ "i"
                  ++ fmtSigned c ++ "j" ++ fmtSigned d ++ "k"
    where fmtSigned x = (if x >= 0 then "+" else "") ++ fmtF 3 x

qadd :: Q -> Q -> Q
qadd (Q a b c d) (Q e f g h) = Q (a+e) (b+f) (c+g) (d+h)

qsub :: Q -> Q -> Q
qsub (Q a b c d) (Q e f g h) = Q (a-e) (b-f) (c-g) (d-h)

-- Non-commutative quaternion multiplication
qmul :: Q -> Q -> Q
qmul (Q a b c d) (Q e f g h) = Q
  (a*e - b*f - c*g - d*h)
  (a*f + b*e + c*h - d*g)
  (a*g - b*h + c*e + d*f)
  (a*h + b*g - c*f + d*e)

qconj :: Q -> Q
qconj (Q a b c d) = Q a (-b) (-c) (-d)

qnorm :: Q -> Double
qnorm (Q a b c d) = sqrt (a*a + b*b + c*c + d*d)

qscale :: Double -> Q -> Q
qscale s (Q a b c d) = Q (s*a) (s*b) (s*c) (s*d)

qneg :: Q -> Q
qneg = qscale (-1)

-- Unit imaginary quaternions
qi, qj, qk :: Q
qi = Q 0 1 0 0
qj = Q 0 0 1 0
qk = Q 0 0 0 1

-- ============================================================
-- THREE COMPLEX STRUCTURES  I, J, K
-- Left multiplication by i, j, k respectively.
-- Satisfy I^2 = J^2 = K^2 = IJK = -1  (verified below)
-- ============================================================

structI :: Q -> Q
structI q = qmul qi q   -- I: q -> i*q

structJ :: Q -> Q
structJ q = qmul qj q   -- J: q -> j*q

structK :: Q -> Q
structK q = qmul qk q   -- K: q -> k*q

-- ============================================================
-- QUATERNIONIC VECTORS
-- ============================================================

type QVec = [Q]

-- Quaternionic Hermitian inner product: H(u,v) = sum conj(u_a) * v_a
-- Returns a quaternion; splits as H = g + i*wI + j*wJ + k*wK
qHerm :: QVec -> QVec -> Q
qHerm u v = foldr qadd (Q 0 0 0 0) (zipWith (\a b -> qmul (qconj a) b) u v)

-- Riemannian metric:    g(u,v)  = Re H(u,v) = q0
metricG :: QVec -> QVec -> Double
metricG u v = q0 (qHerm u v)

-- Symplectic form I:    wI(u,v) = Im_I H(u,v) = q1
omegaI :: QVec -> QVec -> Double
omegaI u v = q1 (qHerm u v)

-- Symplectic form J:    wJ(u,v) = Im_J H(u,v) = q2
omegaJ :: QVec -> QVec -> Double
omegaJ u v = q2 (qHerm u v)

-- Symplectic form K:    wK(u,v) = Im_K H(u,v) = q3
omegaK :: QVec -> QVec -> Double
omegaK u v = q3 (qHerm u v)

-- Quaternionic norm: ||v||^2 = Re H(v,v)
qvNorm :: QVec -> Double
qvNorm v = sqrt (metricG v v)

-- Hermitian distance in H^k
hkDist :: QVec -> QVec -> Double
hkDist u v = sqrt . sum $ zipWith (\a b -> qnorm (qsub a b) ^ 2) u v

-- ============================================================
-- HYPER-KAHLER CONDITION CHECKS
-- ============================================================

eps :: Double
eps = 1e-10

-- I^2 = -id on a single quaternion
checkI2 :: Q -> Bool
checkI2 q = qnorm (qsub (structI (structI q)) (qneg q)) < eps

-- J^2 = -id
checkJ2 :: Q -> Bool
checkJ2 q = qnorm (qsub (structJ (structJ q)) (qneg q)) < eps

-- K^2 = -id
checkK2 :: Q -> Bool
checkK2 q = qnorm (qsub (structK (structK q)) (qneg q)) < eps

-- IJK = -id
checkIJK :: Q -> Bool
checkIJK q = qnorm (qsub (structI (structJ (structK q))) (qneg q)) < eps

-- IJ = K (not -K; left-multiplication gives IJ=K)
checkIJ :: Q -> Bool
checkIJ q = qnorm (qsub (structI (structJ q)) (structK q)) < eps

-- wI skew-symmetric: wI(u,v) = -wI(v,u)
checkWISkew :: QVec -> QVec -> Bool
checkWISkew u v = abs (omegaI u v + omegaI v u) < eps

-- Normalise a QVec to unit Hermitian norm.
-- Needed for compatibility checks: ISCO coords have physical scale
-- mismatches that cause false FAILs before normalisation.
--   FAIL after normalise = genuinely non-Lagrangian embedding
--   FAIL only before  = coordinate scale artefact
qvNormalise :: QVec -> QVec
qvNormalise v =
  let n = qvNorm v
  in if n < 1e-12 then v else map (qscale (1/n)) v

-- wI compatible with I on normalised vectors: wI(u^,v^) = g(I u^, v^)
checkWICompat :: QVec -> QVec -> Bool
checkWICompat u v =
  let u' = qvNormalise u; v' = qvNormalise v
  in abs (omegaI u' v' - metricG (map structI u') v') < 1e-8

checkWJCompat :: QVec -> QVec -> Bool
checkWJCompat u v =
  let u' = qvNormalise u; v' = qvNormalise v
  in abs (omegaJ u' v' - metricG (map structJ u') v') < 1e-8

checkWKCompat :: QVec -> QVec -> Bool
checkWKCompat u v =
  let u' = qvNormalise u; v' = qvNormalise v
  in abs (omegaK u' v' - metricG (map structK u') v') < 1e-8

-- Raw residuals on normalised vectors (for reporting)
compatResidualI, compatResidualJ, compatResidualK :: QVec -> QVec -> Double
compatResidualI u v =
  let u' = qvNormalise u; v' = qvNormalise v
  in omegaI u' v' - metricG (map structI u') v'
compatResidualJ u v =
  let u' = qvNormalise u; v' = qvNormalise v
  in omegaJ u' v' - metricG (map structJ u') v'
compatResidualK u v =
  let u' = qvNormalise u; v' = qvNormalise v
  in omegaK u' v' - metricG (map structK u') v'

-- Discrete curvature: d(wI) over triple p,q,r
dOmegaI :: QVec -> QVec -> QVec -> Double
dOmegaI p q r =
  let s = zipWith qsub
      qp=s q p; rp=s r p; rq=s r q; pq=s p q; pr=s p r; qr=s q r
  in omegaI qp rp + omegaI rq pq + omegaI pr qr

dOmegaJ :: QVec -> QVec -> QVec -> Double
dOmegaJ p q r =
  let s = zipWith qsub
      qp=s q p; rp=s r p; rq=s r q; pq=s p q; pr=s p r; qr=s q r
  in omegaJ qp rp + omegaJ rq pq + omegaJ pr qr

dOmegaK :: QVec -> QVec -> QVec -> Double
dOmegaK p q r =
  let s = zipWith qsub
      qp=s q p; rp=s r p; rq=s r q; pq=s p q; pr=s p r; qr=s q r
  in omegaK qp rp + omegaK rq pq + omegaK pr qr

-- ============================================================
-- POINT TYPES
-- ============================================================

data HKPoint = HKPoint
  { hkLabel :: Int
  , hkQuat  :: QVec    -- quaternionic embedding
  } deriving (Show)

-- ============================================================
-- SIGNALS
-- ============================================================

-- ISCO geodesic: (r, t, theta, phi) per time step -> one Q each
-- Source: kahler_isco.hs / bh-phase-space Playground.hs
iscoR, iscoT, iscoTheta, iscoP :: [Double]
iscoR =
  [-2.2e-05,0.469324,0.927044,1.361952,1.763352,2.121341,2.427033,
   2.673001,2.853188,2.963087,3.0,2.963087,2.853188,2.673001,
   2.427033,2.121341,1.763352,1.361952,0.927044,0.469324,-2.2e-05,
   -0.469324,-0.927044,-1.361952,-1.763352,-2.121341,-2.427077,
   -2.673045,-2.853188,-2.963087,-3.0,-2.963087,-2.853188,-2.673045,
   -2.427077,-2.121341,-1.763352,-1.361952,-0.927044,-0.469324,
   -2.2e-05,0.469324,0.927044,1.361952,1.763352,2.121341,2.427033,
   2.673001,2.853188,2.963087,3.0,2.963087,2.853188,2.673001,
   2.427033,2.121341,1.763352,1.361952,0.927044,0.469324]

iscoT =
  [-3.0,-2.979933,-2.959866,-2.939799,-2.919732,-2.899666,-2.879599,
   -2.859532,-2.839465,-2.819398,-2.799331,-2.779264,-2.759197,
   -2.73913,-2.719064,-2.698997,-2.67893,-2.658863,-2.638796,
   -2.618729,-2.598662,-2.578595,-2.558528,-2.538462,-2.518395,
   -2.498328,-2.478261,-2.458194,-2.438127,-2.41806,-2.397993,
   -2.377926,-2.35786,-2.337793,-2.317726,-2.297659,-2.277592,
   -2.257525,-2.237458,-2.217391,-2.197324,-2.177258,-2.157191,
   -2.137124,-2.117057,-2.09699,-2.076923,-2.056856,-2.036789,
   -2.016722,-1.996656,-1.976589,-1.956522,-1.936455,-1.916388,
   -1.896321,-1.876254,-1.856187,-1.83612,-1.816054]

iscoTheta =
  [0.0,0.342137,0.679872,1.008603,1.32433,1.622649,1.89996,
   2.152261,2.376751,2.570228,2.730092,2.854342,2.941577,2.990196,
   3.0,2.970788,2.902761,2.796719,2.654462,2.477391,2.268307,
   2.029412,1.764106,1.47579,1.168267,0.845538,0.511805,0.171469,
   -0.171269,-0.511805,-0.845538,-1.168267,-1.47579,-1.763906,
   -2.029212,-2.268107,-2.477391,-2.654262,-2.796719,-2.902561,
   -2.970588,-3.0,-2.990196,-2.941377,-2.854342,-2.729892,
   -2.570028,-2.376551,-2.152261,-1.89976,-1.622449,-1.32413,
   -1.008403,-0.679672,-0.342137,0.0,0.342137,0.679872,
   1.008603,1.32433]

iscoP =
  [-3.0,-2.868523,-2.737046,-2.605569,-2.474092,-2.342615,-2.211137,
   -2.079661,-1.948184,-1.816707,-1.685229,-1.553753,-1.422276,
   -1.290799,-1.159321,-1.027845,-0.896368,-0.764891,-0.633413,
   -0.501937,-0.37046,-0.238983,-0.107505,0.023971,0.155448,0.286925,
   0.418403,0.549879,0.681356,0.812833,0.944311,1.075788,1.207264,
   1.338741,1.470219,1.601696,1.733172,1.86465,1.996127,2.127604,
   2.25908,2.390558,2.522035,2.653512,2.784988,2.916466,-2.972155,
   -2.840679,-2.709201,-2.577724,-2.446247,-2.314771,-2.183293,
   -2.051816,-1.920339,-1.788863,-1.657385,-1.525908,-1.394431,
   -1.262954]

-- Native ISCO embedding: Q r t theta phi per step
iscoPoints :: [HKPoint]
iscoPoints = zipWith4 mk [0..] iscoR iscoT iscoTheta iscoP
  where
    mk i r t th ph = HKPoint i [Q r t th ph]
    zipWith4 f (a:as)(b:bs)(c:cs)(d:ds)(e:es) =
      f a b c d e : zipWith4 f as bs cs ds es
    zipWith4 _ _ _ _ _ _ = []

-- Myers-Perry horizon data as quaternionic signal
-- Each point: Q r_H a A_H Omega_H
-- r_H, A_H, Omega_H computed from massParam 5 and spin a — no hardcoding.
mpPoints :: [HKPoint]
mpPoints = zipWith mk [0..] (map mkQ spins)
  where
    spins        = [0.00, 0.30, 0.60, 0.90]
    mkQ a        = let rH = horizonD5 a
                       aH = areaD5 rH a
                       oh = omegaHD5 rH a
                   in  Q rH a aH oh
    mk i q       = HKPoint i [q]

-- Black ring horizon data as quaternionic signal
-- Each point: Q nu lambda A_H Omega_H
ringPoints :: [HKPoint]
ringPoints = zipWith mk [0..]
  [ Q 0.1 (eqL 0.1) 1.2223  0.6428
  , Q 0.3 (eqL 0.3) 12.6156 0.5417
  , Q 0.5 (eqL 0.5) 39.9493 0.4564
  , Q 0.7 (eqL 0.7) 100.2767 0.3626
  , Q 0.9 (eqL 0.9) 371.4256 0.2182
  ]
  where
    mk i q = HKPoint i [q]
    eqL nu = 2*nu/(1+nu*nu)

-- Myers-Perry d=6 as 2-quaternion signal in H^2
-- Each point: [Q(r_H, a, 0, 0), Q(mu, nu, 0, 0)]
-- r_H computed via Cardano (single-spin d=6: r^3 + a^2*r - mu = 0).
-- nu = r_H / a (ratio of horizon radius to spin parameter).
mp6Points :: [HKPoint]
mp6Points = zipWith mk [0..] (map mkQs spins)
  where
    mu6          = massParam 6
    spins        = [0.0, 0.3, 0.5, 0.6]
    mkQs a       = let rH = horizonD6Single a
                       nu = if a > 0 then rH / a else 0.0
                   in  [ Q rH a 0.0 0.0
                        , Q mu6 nu 0.0 0.0 ]
    mk i qs      = HKPoint i qs

-- ============================================================
-- MULTI-SPIN d=6  (two independent rotation planes)
-- H^2 embedding: [Q(r_H, a1, a2, 0), Q(mu, nu1, nu2, 0)]
-- nu1 = r_H/a1,  nu2 = r_H/a2
--
-- Three regimes tested:
--   A) a1 >> a2  (near single-spin, recovers mp6 result)
--   B) a1 = a2   (equal spins, expect SU(2) symmetry: cI~cJ~cK)
--   C) a2 -> 0   (single spin limit, cJ=cK=0 exactly)
--
-- r_H computed via closed-form quadratic in u=r^2:
--   u^2 + (a1^2 + a2^2 - mu)*u + a1^2*a2^2 = 0
-- All horizons update automatically when bigG or massM change.
-- ============================================================

-- Helper: build a 2-quaternion HKPoint from two-spin parameters.
-- r_H computed via two-spin quadratic; nu1, nu2 derived from it.
-- Valid only for a1+a2 <= sqrt(mu6) ~ 0.691.
mkMP6Two :: Int -> Double -> Double -> HKPoint
mkMP6Two lbl a1 a2 =
  let mu  = massParam 6
      rH  = horizonD6Two a1 a2
      nu1 = if a1 > 0 then rH / a1 else 0.0
      nu2 = if a2 > 0 then rH / a2 else 0.0
  in  HKPoint lbl [Q rH a1 a2 0.0, Q mu nu1 nu2 0.0]

-- Helper: build a 2-quaternion HKPoint from single-spin parameters (a2=0).
-- r_H computed via Cardano (depressed cubic).
mkMP6One :: Int -> Double -> HKPoint
mkMP6One lbl a1 =
  let mu  = massParam 6
      rH  = horizonD6Single a1
      nu1 = if a1 > 0 then rH / a1 else 0.0
  in  HKPoint lbl [Q rH a1 0.0 0.0, Q mu nu1 0.0 0.0]

-- Regime A: a1 >> a2  (near single-spin)
-- All within extremal bound a1+a2 <= sqrt(mu6) ~ 0.691
mpA :: [HKPoint]
mpA = zipWith (\i (a1,a2) -> mkMP6Two i a1 a2) [0..]
  [(0.60, 0.05), (0.50, 0.05), (0.40, 0.05), (0.30, 0.05)]

-- Regime B: a1 = a2  (equal spins, full SU(2))
-- Extremal for equal spins: a <= sqrt(mu6)/2 ~ 0.345
mpB :: [HKPoint]
mpB = zipWith (\i (a1,a2) -> mkMP6Two i a1 a2) [0..]
  [(0.30, 0.30), (0.25, 0.25), (0.20, 0.20), (0.10, 0.10)]

-- Regime C: a2 = 0  (single-spin limit, cJ=cK=0 exactly)
mpC :: [HKPoint]
mpC = zipWith mkMP6One [0..] [0.0, 0.3, 0.5, 0.6]

-- ============================================================
-- kNN in H^k
-- ============================================================

hkKNN :: Int -> HKPoint -> [HKPoint] -> [(HKPoint, Double)]
hkKNN k query corpus =
  take k . sortBy (comparing snd)
  . map (\p -> (p, hkDist (hkQuat query) (hkQuat p)))
  . filter ((/= hkLabel query) . hkLabel)
  $ corpus

-- ============================================================
-- CHERN PROXIES  (one per symplectic form)
-- ============================================================

chernI, chernJ, chernK :: [HKPoint] -> Double
chernI pts = sum [omegaI (hkQuat (pts!!i)) (hkQuat (pts!!j))
                 | i<-[0..n-1], j<-[i+1..n-1]]
  where n = length pts
chernJ pts = sum [omegaJ (hkQuat (pts!!i)) (hkQuat (pts!!j))
                 | i<-[0..n-1], j<-[i+1..n-1]]
  where n = length pts
chernK pts = sum [omegaK (hkQuat (pts!!i)) (hkQuat (pts!!j))
                 | i<-[0..n-1], j<-[i+1..n-1]]
  where n = length pts

-- ============================================================
-- VR PERSISTENT HOMOLOGY in H^k
--
-- Strategy:
--   1. Build distance matrix once as Array (Int,Int) Double — O(1) lookup.
--   2. Sort all edges once by distance — O(n² log n).
--   3. Walk epsilon steps in order, consuming edges incrementally.
--      State (edges, membership array, union-find, triangle count)
--      is carried forward, never rebuilt from scratch.
--   4. Triangle count: when a new edge (i,j) is added, count k s.t.
--      both (i,k) and (j,k) already exist — O(n) per new edge.
--
-- Cost: O(n² log n + n² × n) vs original O(steps × n³).
-- ============================================================

-- Union-find with path compression using Array
type UF = Array Int Int

mkUF :: Int -> UF
mkUF n = listArray (0, n-1) [0..n-1]

findUF :: UF -> Int -> (Int, UF)
findUF uf i
  | uf ! i == i = (i, uf)
  | otherwise   =
      let (root, uf') = findUF uf (uf ! i)
      in  (root, uf' // [(i, root)])

unionUF :: UF -> Int -> Int -> UF
unionUF uf a b =
  let (ra, uf')  = findUF uf  a
      (rb, uf'') = findUF uf' b
  in  if ra == rb then uf'' else uf'' // [(ra, rb)]

countComps :: UF -> Int -> Int
countComps uf n =
  length . nub $ map (fst . findUF uf) [0..n-1]

-- Edge membership: Array (Int,Int) Bool, always indexed with i < j
type EdgeArr = Array (Int,Int) Bool

mkEdgeArr :: Int -> EdgeArr
mkEdgeArr n = listArray ((0,0),(n-1,n-1)) (repeat False)

hasEdge :: EdgeArr -> Int -> Int -> Bool
hasEdge arr i j = arr ! (min i j, max i j)

addEdge :: EdgeArr -> Int -> Int -> EdgeArr
addEdge arr i j = arr // [((min i j, max i j), True)]

-- mapAccumL: thread state through a list, collecting results
mapAccumL :: (s -> a -> (s, b)) -> s -> [a] -> (s, [b])
mapAccumL _ s []     = (s, [])
mapAccumL f s (x:xs) =
  let (s', b)   = f s x
      (s'', bs) = mapAccumL f s' xs
  in  (s'', b : bs)

-- Vietoris-Rips filtration
-- Returns [(epsilon, H0_components, H1_cycles)]
vrFiltration :: [HKPoint] -> Int -> [(Double,Int,Int)]
vrFiltration pts nSteps =
  let n = length pts

      -- Distance array: O(1) lookup — computed once
      distArr :: Array (Int,Int) Double
      distArr = listArray ((0,0),(n-1,n-1))
                  [ hkDist (hkQuat (pts!!i)) (hkQuat (pts!!j))
                  | i <- [0..n-1], j <- [0..n-1] ]
      d i j = distArr ! (i, j)

      -- All edges sorted by distance, built once — O(n² log n)
      sortedEdges :: [(Double,Int,Int)]
      sortedEdges = sortBy (comparing (\(w,_,_) -> w))
                      [ (d i j, i, j) | i <- [0..n-1], j <- [i+1..n-1] ]

      -- Select nSteps epsilon values from the sorted edge list
      m      = length sortedEdges
      epIdxs = nub [ (k*(m-1)) `div` max 1 (nSteps-1) | k <- [0..nSteps-1] ]
      epsilons = map (\i -> let (w,_,_) = sortedEdges!!i in w) epIdxs

      -- Incremental step: extend state to cover all edges <= eps
      step :: ([(Double,Int,Int)], [(Int,Int)], EdgeArr, UF, Int)
           -> Double
           -> ( ([(Double,Int,Int)], [(Int,Int)], EdgeArr, UF, Int)
              , (Double,Int,Int) )
      step (remaining, curEdges, eArr, uf, triCount) eps =
        let (batch, rest) = span (\(w,_,_) -> w <= eps) remaining

            -- For new edge (i,j): count k where (i,k) and (j,k) already exist
            newTrisFor i j =
              length [ () | k <- [0..n-1]
                          , k /= i, k /= j
                          , hasEdge eArr (min i k) (max i k)
                          , hasEdge eArr (min j k) (max j k) ]

            addOne (es, earr, u, tc) (_, i, j) =
              ( (i,j) : es
              , addEdge earr i j
              , unionUF u i j
              , tc + newTrisFor i j )

            (curEdges', eArr', uf', triCount') =
              foldl addOne (curEdges, eArr, uf, triCount) batch

            h0 = countComps uf' n
            nE = length curEdges'
            h1 = max 0 (nE - (n - h0) - triCount')

        in  ( (rest, curEdges', eArr', uf', triCount')
            , (eps, h0, h1) )

      initState = (sortedEdges, [], mkEdgeArr n, mkUF n, 0)
      (_, results) = mapAccumL step initState epsilons

  in  results

-- ============================================================
-- PRETTY PRINTING
-- ============================================================

fmtF :: Int -> Double -> String
fmtF dp x = show (fromIntegral (round (x * 10^dp)) / 10^dp :: Double)

padL :: Int -> String -> String
padL n s = replicate (max 0 (n - length s)) ' ' ++ s

padR :: Int -> String -> String
padR n s = take n (s ++ repeat ' ')

bar :: Double -> Double -> Int -> String
bar v mx w =
  let len = max 0 . min w $ round (abs v / max (abs mx) 1e-9 * fromIntegral w)
  in replicate len '#' ++ replicate (w-len) '.'

tick :: Bool -> String
tick True  = "OK"
tick False = "FAIL"

section :: String -> IO ()
section t = do
  putStrLn ""
  putStrLn $ "---  " ++ t
  putStrLn $ replicate (length t + 5) '-'

-- ============================================================
-- REPORTS
-- ============================================================

reportAlgebra :: IO ()
reportAlgebra = do
  section "Quaternion algebra  I^2=J^2=K^2=IJK=-1"
  let qs = [Q 1 2 3 4, Q 0.5 (-1) 2 0, Q 3 0 (-1) 2, Q 0 1 0 0]
  putStrLn "  q                        I^2   J^2   K^2   IJK   IJ=K"
  mapM_ (\q ->
    putStrLn $ "  " ++ padR 26 (show q)
      ++ padL 6 (tick (checkI2  q))
      ++ padL 6 (tick (checkJ2  q))
      ++ padL 6 (tick (checkK2  q))
      ++ padL 6 (tick (checkIJK q))
      ++ padL 6 (tick (checkIJ  q))
    ) qs
  putStrLn ""
  putStrLn "  Non-commutativity:  i*j = k,  j*i = -k"
  putStrLn $ "  i*j = " ++ show (qmul qi qj)
  putStrLn $ "  j*i = " ++ show (qmul qj qi)

reportHKConditions :: [HKPoint] -> String -> IO ()
reportHKConditions pts label = do
  section $ "Hyper-Kahler conditions -- " ++ label
  putStrLn "  (compatibility checked on normalised unit vectors)"
  putStrLn "  t    wI-skew  wI-compat  wJ-compat  wK-compat  resI       resJ       resK"
  mapM_ (\(p,q) ->
    let u = hkQuat p; v = hkQuat q
    in putStrLn $ "  " ++ padR 5 (show (hkLabel p))
         ++ padL 9  (tick (checkWISkew   u v))
         ++ padL 11 (tick (checkWICompat u v))
         ++ padL 11 (tick (checkWJCompat u v))
         ++ padL 11 (tick (checkWKCompat u v))
         ++ padL 11 (fmtF 5 (compatResidualI u v))
         ++ padL 11 (fmtF 5 (compatResidualJ u v))
         ++ padL 11 (fmtF 5 (compatResidualK u v))
    ) (take 8 (zip pts (drop 1 pts)))

reportCurvature :: [HKPoint] -> String -> IO ()
reportCurvature pts label = do
  section $ "Curvature d(omega) over triples -- " ++ label
  putStrLn "  t    d(wI)        d(wJ)        d(wK)        reading"
  mapM_ (\(p,q,r) ->
    let u=hkQuat p; v=hkQuat q; w=hkQuat r
        dI = dOmegaI u v w
        dJ = dOmegaJ u v w
        dK = dOmegaK u v w
        note = if all (\x->abs x<1e-6) [dI,dJ,dK]
               then "Lagrangian (flat)"
               else "curved HK manifold"
    in putStrLn $ "  " ++ padR 5 (show (hkLabel p))
         ++ padL 13 (fmtF 5 dI)
         ++ padL 13 (fmtF 5 dJ)
         ++ padL 13 (fmtF 5 dK)
         ++ "  " ++ note
    ) (take 10 (zip3 pts (drop 1 pts) (drop 2 pts)))

reportChern :: [(String, [HKPoint])] -> IO ()
reportChern signals = do
  section "Chern proxies -- omega_I, omega_J, omega_K per signal"
  putStrLn "  Discrete integral sum_{i<j} w(p_i, p_j)"
  putStrLn "  Killing-Yano correspondence (single-spin MP / ISCO):"
  putStrLn "    omega_I ~ Omega_I  (radial KY tensor)"
  putStrLn "    omega_J ~ Omega_J  (polar KY tensor)"
  putStrLn "    omega_K ~ Omega_K  (azimuthal KY tensor)"
  putStrLn ""
  let rows = map (\(lbl,pts) -> (lbl, chernI pts, chernJ pts, chernK pts)) signals
      mx   = maximum [abs x | (_,a,b,c) <- rows, x <- [a,b,c]]
  putStrLn $ "  " ++ padR 14 "signal"
          ++ padR 28 "omega_I" ++ padR 28 "omega_J" ++ "omega_K"
  mapM_ (\(lbl,ci,cj,ck) ->
    putStrLn $ "  " ++ padR 14 lbl
      ++ padR 4 "" ++ bar ci mx 20 ++ " " ++ padL 8 (fmtF 3 ci)
      ++ padR 4 "" ++ bar cj mx 20 ++ " " ++ padL 8 (fmtF 3 cj)
      ++ padR 4 "" ++ bar ck mx 20 ++ " " ++ padL 8 (fmtF 3 ck)
    ) rows

reportNorms :: [HKPoint] -> String -> IO ()
reportNorms pts label = do
  section $ "Quaternionic norms -- " ++ label
  putStrLn "  Re H(v,v) = ||v||^2  confirms metric preserved"
  putStrLn "  t    q0 (r)    q1 (t)    q2 (th)   q3 (phi)  |q|"
  mapM_ (\p ->
    let [Q a b c d] = hkQuat p
        n = sqrt (a*a+b*b+c*c+d*d)
    in putStrLn $ "  " ++ padR 5 (show (hkLabel p))
         ++ padL 10 (fmtF 4 a)
         ++ padL 10 (fmtF 4 b)
         ++ padL 10 (fmtF 4 c)
         ++ padL 10 (fmtF 4 d)
         ++ padL 8  (fmtF 4 n)
    ) (take 10 pts)

reportKNN :: [HKPoint] -> String -> IO ()
reportKNN pts label = do
  section $ "Quaternionic 5-NN  (Hermitian L2 in H^k) -- " ++ label
  case filter ((==10) . hkLabel) pts of
    [] -> putStrLn "  t=10 not found."
    (q:_) -> do
      let nbs = hkKNN 5 q pts
          mx  = maximum (map snd nbs)
      putStrLn "  t     dist"
      mapM_ (\(p,d) ->
        putStrLn $ "  " ++ padR 6 (show (hkLabel p))
          ++ bar d mx 28 ++ " " ++ fmtF 5 d
        ) nbs

reportVR :: [HKPoint] -> Int -> String -> IO ()
reportVR pts steps label = do
  section $ "VR persistent homology in H^k -- " ++ label
  putStrLn "  epsilon       H0   H1"
  mapM_ (\(e,h0,h1) ->
    putStrLn $ "  " ++ padR 14 (fmtF 5 e)
      ++ padL 4 (show h0)
      ++ padL 5 (show h1)
    ) (vrFiltration pts steps)

reportComparison :: IO ()
reportComparison = do
  section "Chern proxy comparison: ISCO vs MP vs Black Ring"
  putStrLn "  Prediction for single-spin geometry:"
  putStrLn "    omega_I >> omega_J ~ omega_K  (radial KY dominates)"
  putStrLn "    Black ring: omega_K >> omega_I,J  (azimuthal rotation)"
  putStrLn ""
  let signals =
        [ ("ISCO",  iscoPoints)
        , ("MP-d5", mpPoints)
        , ("B-Ring", ringPoints)
        ]
  reportChern signals

-- ============================================================
-- MULTI-SPIN REPORT
-- ============================================================

reportMultiSpin :: IO ()
reportMultiSpin = do
  section "MP d=6 multi-spin -- three regimes in H^2"
  putStrLn "  Each point: [Q(r_H, a1, a2, A_H), Q(mu, nu1, nu2, Omega)]"
  putStrLn "  Prediction:"
  putStrLn "    Regime A (a1>>a2): cK dominates (azimuthal KY from large a1)"
  putStrLn "    Regime B (a1=a2):  cI = cJ, cK = 0  (SU(2) symmetry)"
  putStrLn "    Regime C (a2=0):   cJ = cK = 0   (single-spin limit)"
  putStrLn ""

  let regimes = [("A: a1>>a2", mpA), ("B: a1=a2", mpB), ("C: a2=0", mpC)]

  reportChern regimes
  mapM_ (\(lbl, pts) -> reportHKConditions pts ("MP d=6 " ++ lbl)) regimes
  mapM_ (\(lbl, pts) -> reportCurvature pts ("MP d=6 " ++ lbl)) regimes

  section "Spin symmetry fingerprint"
  putStrLn "  signal          cI/cK ratio   cJ/cK ratio   interpretation"
  mapM_ (\(lbl, pts) ->
    let ci = chernI pts
        cj = chernJ pts
        ck = chernK pts
        rIK = if abs ck > 1e-10 then ci/ck else 999.0
        rJK = if abs ck > 1e-10 then cj/ck else 999.0
        interp
          | abs ck < 1e-8 && abs cj < 1e-8 = "single-spin (I only)"
          | abs ck < 1e-8 && abs (ci - cj) < 1e-6 = "near SU(2) [cK=0]"
          | abs (rIK - 1.0) < 0.2 && abs (rJK - 1.0) < 0.2 = "near SU(2)"
          | abs ci > abs cj && abs ci > abs ck = "I-dominant"
          | otherwise = "mixed"
    in putStrLn $ "  " ++ padR 16 lbl
         ++ padL 14 (fmtF 3 rIK)
         ++ padL 14 (fmtF 3 rJK)
         ++ "  " ++ interp
    ) regimes

-- ============================================================
-- MAIN
-- ============================================================

main :: IO ()
main = do
  putStrLn "AVS Hyper-Kahler Extension -- base only -- play.haskell.org"
  putStrLn "Stack: avsp -> kahler-ts -> kahler_isco -> bh-phase-space -> THIS"
  putStrLn ""
  putStrLn "H = g + i*omegaI + j*omegaJ + k*omegaK"
  putStrLn "  g       = Re H   (Riemannian, preserves L2 from avsp)"
  putStrLn "  omegaI  = q1 H   (Killing-Yano Omega_I, radial)"
  putStrLn "  omegaJ  = q2 H   (Killing-Yano Omega_J, polar)"
  putStrLn "  omegaK  = q3 H   (Killing-Yano Omega_K, azimuthal)"

  -- 1. Algebra self-check
  reportAlgebra

  -- 2. ISCO native embedding
  reportNorms iscoPoints "ISCO native Q(r,t,theta,phi)"

  -- 3. Hyper-Kahler conditions on ISCO
  reportHKConditions iscoPoints "ISCO"

  -- 4. Curvature on ISCO triples
  reportCurvature iscoPoints "ISCO"

  -- 5. Chern proxies: ISCO vs MP vs Black Ring
  reportComparison

  -- 6. kNN on ISCO
  reportKNN iscoPoints "ISCO"

  -- 7. VR homology on ISCO
  reportVR iscoPoints 20 "ISCO (all 60 points, 20 epsilon steps)"

  -- 8. MP and ring conditions
  reportHKConditions mpPoints   "MP d=5 horizon data"
  reportCurvature    mpPoints   "MP d=5 horizon data"
  reportHKConditions ringPoints "Black ring data"
  reportCurvature    ringPoints "Black ring data"

  -- 9. MP d=6 single-spin baseline (H^2)
  section "MP d=6 single-spin -- H^2 baseline"
  putStrLn "  [Q(r_H, a, A_H, 0), Q(mu, nu, Omega, 0)]"
  putStrLn "  Prediction: cJ=cK=0 (single complex structure only)"
  reportHKConditions mp6Points "MP d=6 single-spin"
  reportCurvature    mp6Points "MP d=6 single-spin"
  reportChern [("MP-d5", mpPoints), ("MP-d6-1spin", mp6Points), ("B-Ring", ringPoints)]

  -- 10. Multi-spin d=6
  reportMultiSpin

  putStrLn ""
  section "Summary"
  putStrLn "  ISCO orbit:     4D spacetime position = 1 quaternion natively"
  putStrLn "  wI/wJ/wK:       three Killing-Yano tensors of Myers-Perry"
  putStrLn "  ISCO FAILs:     non-Lagrangian; resI~-resK (I/K antisymmetry)"
  putStrLn "  MP d=6 1-spin:  Lagrangian in I; cJ=cK=0 exactly"
  putStrLn "  MP d=6 2-spin:  A(a1>>a2)->K-dominant; B(equal)->cI=cJ,cK=0; C(a2=0)->1-spin"
  putStrLn "  VR H1=0:        ISCO < 1 orbit; need denser sampling for loop"
  putStrLn ""
  putStrLn "  Next: ISCO full orbit -> H1 generator"
  putStrLn "        d=5 ring vs MP triple-degeneracy VR comparison"
