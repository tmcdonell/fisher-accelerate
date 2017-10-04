{-# LANGUAGE ConstraintKinds   #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE ViewPatterns      #-}

-- This module will only define functions in Accelerate, so sometimes it is
-- useful to use the following two extensions. In particular, this allows us to
-- use the usual if-then-else syntax with Accelerate Acc and Exp terms, rather
-- than `(|?)` and `(?)` respectively.
--
{-# LANGUAGE RebindableSyntax  #-}
{-# LANGUAGE NoImplicitPrelude #-}

module Solve where

import Data.Array.Accelerate
import Data.Array.Accelerate.Control.Lens


-- Useful type synonyms
type R        = Double
type Field e  = Array DIM2 e


-- | Take a single time step
--
step :: Acc (Scalar R)
     -> Acc (Field R)
     -> Acc (Scalar Bool, Field R)
step dt x0 =
  let
      -- maximum number of iterations for convergence
      maxiters  = 50
      yes       = unit (constant True)
      no        = unit (constant False)

      -- convergence tolerance
      tolerance = 1.0e-3

      -- loop condition
      cont :: Acc (Scalar Bool, Scalar Int, Field R)
           -> Acc (Scalar Bool)
      cont st = zipWith (\converged i -> not converged && i < maxiters) (st ^. _1) (st ^. _2)

      -- loop body
      body :: Acc (Scalar Bool, Scalar Int, Field R)
           -> Acc (Scalar Bool, Scalar Int, Field R)
      body st =
        let
            it        = st ^. _2
            x_new     = st ^. _3

            -- compute residual
            b         = diffusion dt x0 x_new
            residual  = norm2 (flatten b)

        -- check for convergence
        in
        if the residual < tolerance
          then lift (yes, it, x_new)
          else
            let -- solve linear system to get -delta_x
                (ok, dx)  = unlift $ solve dt x0 b
                x_new'    = zipWith (-) x_new dx
            in
            if the ok
              then lift (no, map (+1) it,   x_new')
              else lift (no, unit maxiters, x_new')

      -- solver loop
      result  = awhile cont body (lift (no, unit 0, x0))
  in
  lift (result ^. _1, result ^. _3)


-- | Conjugate gradient solver routine.
--
-- Solve the linear system \( A * x = b \) for @x@.
--
-- The matrix A is implicit in the objective function for the diffusion
-- equation. The input @x@ is used as the initial guess at the solution.
--
solve :: Acc (Scalar R)
      -> Acc (Field R)
      -> Acc (Field R)
      -> Acc (Scalar Bool, Field R)
solve dt x0 b0 =
  let
      -- epsilon value used for matrix-vector approximation
      eps       = 1.0e-8
      eps_inv   = 1.0 / eps

      -- maximum number of iterations for convergence
      maxiters  = 200

      -- convergence tolerance
      tolerance = 1.0e-3

      -- matrix-vector multiplication is approximated with:
      --
      --   A*v = 1/epsilon * ( F( x+epsilon*v ) - F(x) )
      --       = 1/epsilon * ( F( x+epsilon*v ) - Fx_old )
      --
      fx0       = diffusion dt x0 x0
      v0        = map ((1+eps) *) x0
      fx'       = diffusion dt x0 v0
      r0        = zipWith3 (\b l r -> b - eps_inv * (l-r)) b0 fx' fx0
      rnew0     = let r0_ = flatten r0 in dot r0_ r0_
      success   = map (\u -> sqrt u < tolerance) rnew0

      -- loop condition
      -- continues the loop body while this evaluates to True
      cont :: Acc (Scalar Bool, Scalar Int, Field R, Field R, Field R, Scalar R)
           -> Acc (Scalar Bool)
      cont st = zipWith (\converged i -> not converged && i < maxiters) (st ^. _1) (st ^. _2)

      -- loop body
      -- keep executing this until the parameters converge (or the limit is reached)
      body :: Acc (Scalar Bool, Scalar Int, Field R, Field R, Field R, Scalar R)
           -> Acc (Scalar Bool, Scalar Int, Field R, Field R, Field R, Scalar R)
      body st =
        let
            -- extract loop-carried variables
            --
            -- We can also use 'unlift' to convert the Acc-of-tuple into
            -- a tuple-of-Acc (i.e. Acc (a,b,c) -> (Acc a, Acc b, Acc c)).
            --
            -- Note however that 'unlift' can be difficult to type-check unless
            -- you use every part of the result. To avoid this, we generally
            -- must give a type signature to the result of 'unlift'.
            -- Alternatively, and as we have done here, we can use operators
            -- from the 'lens-accelerate' package to peek at individual
            -- components of the tuple, which avoid the type inference problem.
            --
            iter      = st ^. _2
            p         = st ^. _3
            r         = st ^. _4
            x         = st ^. _5
            rold      = st ^. _6

            -- Ap = A * p
            v         = zipWith (+) x0 (map (*eps) p)
            fx        = diffusion dt x0 v
            ap        = map (*eps_inv) (zipWith (-) fx fx0)

            alpha     = the rold / the (dot (flatten p) (flatten ap))

            x'        = zipWith (+) x (map (* alpha) p)
            r'        = zipWith (-) r (map (* alpha) ap)

            -- find the new norm
            rnew      = let r'_ = flatten r' in dot r'_ r'_

            -- test for convergence
            converged = map (\u -> sqrt u < tolerance) rnew
            iter'     = map (+1) iter

            -- p is the direction of convergence?
            p'        = if the converged
                          then p
                          else zipWith (+) r' (map (* (the rnew / the rold)) p)
        in
        lift (converged, iter', p', r', x', rnew)

      -- solver loop
      result  = awhile cont body (lift (success, unit 1, r0, r0, x0, rnew0))
  in
  lift (result ^. _1, result ^. _5)


-- | Compute the inner product of two vectors @x@ and @y@.
--
dot :: Num e => Acc (Vector e) -> Acc (Vector e) -> Acc (Scalar e)
dot xs ys
  = fold (+) 0
  $ zipWith (*) xs ys


-- | Compute the L^2 norm (Euclidean distance) of the given vector
--
norm2 :: Floating e => Acc (Vector e) -> Acc (Scalar e)
norm2 xs
  = map sqrt
  $ dot xs xs


-- | The problem is over a rectangular grid of (nx * ny) points. A finite volume
-- discritisation and method of lines gives the following ordinary differential
-- equation for each grid point:
--
-- \[
--    \frac{d s_{ij}}{d t} = \frac{D}{\Delta x^2} [ -4s_{ij} + s_{i-1,j} + s_{i+1,j} + s_{i,j-1} + s_{i,j+1}  ] + R s_{ij}(1 - s_{ij})
-- \]
--
-- Which is expressed as the following nonlinear problem:
--
-- \[
--    f_{ij} = [ -(4 + \alpha)s_{ij} + s_{i-1,j} + s_{i+1,j} + s_{i,j-1} + s_{i,j+1} + \beta s_{ij}(1 - s_{ij}) ]^{k+1} + \alpha s_{ij}^k = 0
-- \]
--
diffusion
    :: Acc (Scalar R)
    -> Acc (Field R)
    -> Acc (Field R)
    -> Acc (Field R)
diffusion (the -> dt) x0 x1 = stencil2 f dirichlet x0 dirichlet x1
  where
    -- Accelerate uses nested tuples to specify stencil patterns, which is
    -- particularly useful for 1D and 2D stencils as they visually represent
    -- which elements of the array are being accessed. Note that we can use as
    -- many or as few of the elements of the stencil as we wish.
    --
    f :: Stencil3x3 R -> Stencil3x3 R -> Exp R
    f st0 st1 = (-(4.0 + alpha) * c) + n + s + e + w
              + alpha * k
              + beta * c * (1 - c)
      where
        ((_,_,_),
         (_,k,_),
         (_,_,_)) = st0

        ((_,n,_),
         (w,c,e),
         (_,s,_)) = st1

    -- For the stencil operator, we need to specify what happens when the
    -- neighbouring points read by the stencil are out-of-bounds.
    --
    -- For this problem we just set all out-of-bounds points to zero. The
    -- function we use here is given out-of-bounds index, which you can examine
    -- to determine which boundary it is. In this way you can, for example,
    -- specify different boundary conditions along different edges, or read in
    -- the boundary values from a separate array.
    --
    -- Other inbuilt boundary conditions include 'wrap', 'mirror', and 'clamp'.
    --
    dirichlet     = function (\_ -> 0)

    -- compute some constants
    Z :. ny :. nx = unlift (shape x0)

    dx            = 1.0 / (fromIntegral (nx - 1))
    dy            = 1.0 / (fromIntegral (ny - 1))
    dxy           = dx * dy
    alpha         = dxy / (1.0 * dt)            -- diffusion coefficient D = 1.0
    beta          = 1000.0 * dxy

