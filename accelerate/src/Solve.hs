{-# LANGUAGE ConstraintKinds   #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE ViewPatterns      #-}

-- This module will only define functions in Accelerate, so sometimes it is
-- useful to use the following two extensions. In particular, this allows us to
-- use the usual if-then-else syntax with Accelerate Acc and Exp terms, rather
-- than `(|?)` and `(?)` (or `acond` and `cond`) respectively.
--
{-# LANGUAGE RebindableSyntax  #-}
{-# LANGUAGE NoImplicitPrelude #-}

module Solve where

import Data.Array.Accelerate
import Data.Array.Accelerate.Control.Lens


-- Useful type synonyms
--
type R        = Double
type Field e  = Array DIM2 e


-- | Take a single time step. Repeatedly call this function until we get to the
-- final solution.
--
-- This corresponds to the main loop of the program in 'c/main.c' line 155.
--
step :: Acc (Scalar R)
     -> Acc (Field R)
     -> Acc (Scalar Bool, Field R)
step dt x0 =
  let
      -- maximum number of iterations for convergence
      maxiters  = 50

      -- convergence tolerance
      tolerance = 1.0e-3

      -- helper constants
      yes       = unit (constant True)
      no        = unit (constant False)

      -- Initial conditions for the loop
      --
      -- This consists of a flag indicating that the solution has converged (and
      -- we can stop), the number of iterations we have taken so far, and the
      -- current solution.
      --
      start :: Acc (Scalar Bool, Scalar Int, Field R)
      start = lift (no, unit 0, x0)

      -- Loop condition
      --
      -- The loop will keep executing while this returns True. We stop once the
      -- solver has converged, or we reach the iteration limit.
      --
      cont :: Acc (Scalar Bool, Scalar Int, Field R)
           -> Acc (Scalar Bool)
      cont state = zipWith (\converged i -> not converged && i < maxiters) (state ^. _1) (state ^. _2)

      -- Loop body
      --
      body :: Acc (Scalar Bool, Scalar Int, Field R)
           -> Acc (Scalar Bool, Scalar Int, Field R)
      body state =
        let
            -- Unpack the loop-carried state variables.
            --
            -- In Accelerate we can use the functions `lift` and `unlift` to
            -- push the `Exp` and `Acc` type through constructors, such as
            -- tuples. That is, we use `unlift` here to convert our Acc-tuple
            -- into a tuple-of-Acc. `lift` is the inverse.
            --
            -- Note that, because we do not use one of the results of the
            -- `unlift`, the type checker will complain that it does not know
            -- what the type of that element should be; thus, we needed to add
            -- an explicit type signature to the result of the `unlift`.
            --
            -- Another approach (and as seen in the `cont` function) is to use
            -- lenses to access each component of the structure. These are
            -- provided by the `lens-accelerate` package. Example:
            --
            -- > it        = state ^. _2
            -- > x_new     = state ^. _3
            --
            (_, it, x_new)  = unlift state  :: (Acc (Scalar Bool), Acc (Scalar Int), Acc (Field R))

            -- Compute residual
            --
            b         = diffusion dt x0 x_new
            residual  = norm2 (flatten b)

        -- Check for convergence
        --
        -- The RebindableSyntax extension lets us use the regular if-then-else
        -- syntax for Accelerate terms, at both the `Acc` and `Exp` level. This
        -- gets converted into the regular function calls `acond` and `cond`
        -- respectively (there are also infix versions of these functions;
        -- `(?|)` and `(?)`). Check the type of those operators to see what the
        -- if-then-else block expects.
        --
        -- Notice that the branches use the `lift` function to tuple-up their
        -- results; this converts the tuple-of-Acc into an Acc-tuple. In this
        -- instance, we have that;
        --
        -- > lift :: (Acc (Scalar Bool), Acc (Scalar Int), Acc (Field R)) -> Acc (Scalar Bool, Scalar Int, Field R)
        --
        in
        if the residual < tolerance
          then lift (yes, it, x_new)                -- solution converged; exit with success
          else
            let -- Solve linear system to get -delta_x
                --
                -- This also tells us (via `ok`) whether the solver converged.
                -- If it did not, then we should also exit this loop to indicate
                -- that the solution has failed.
                --
                -- Notice that we didn't have to add a type signature to
                -- `unlift` here because we use every component of the result,
                -- so, GHC can figure out what the types should be.
                --
                (ok, dx)  = unlift $ solve dt x0 b
                x_new'    = zipWith (-) x_new dx
            in
            if the ok
              then lift (no, map (+1) it,   x_new') -- keep iterating to solution
              else lift (no, unit maxiters, x_new') -- break out of the loop

      -- solver loop
      result  = awhile cont body start
  in
  lift (result ^. _1, result ^. _3)


-- | Conjugate gradient solver routine.
--
-- Solve the linear system \( A * x = b \) for @x@.
--
-- The matrix A is implicit in the objective function for the diffusion
-- equation. The input @x@ is used as the initial guess at the solution.
--
-- This corresponds to the 'ss_cg' function in 'c/linalg.c' line 169.
--
solve :: Acc (Scalar R)
      -> Acc (Field R)
      -> Acc (Field R)
      -> Acc (Scalar Bool, Field R)
solve dt x0 b0 =
  let
      -- epsilon value used for matrix-vector approximation
      eps         = 1.0e-8
      eps_inv     = 1.0 / eps

      -- maximum number of iterations for convergence
      maxiters    = 200

      -- convergence tolerance
      tolerance   = 1.0e-3

      -- matrix-vector multiplication is approximated with:
      --
      --   A*v = 1/epsilon * ( F( x+epsilon*v ) - F(x) )
      --       = 1/epsilon * ( F( x+epsilon*v ) - Fx_old )
      --
      fx_old      = diffusion dt x0 x0
      v0          = map ((1+eps) *) x0
      fx'         = diffusion dt x0 v0
      r0          = zipWith3 (\b l r -> b - eps_inv * (l-r)) b0 fx' fx_old
      rnew0       = let r0_ = flatten r0 in dot r0_ r0_
      converged0  = map (\u -> sqrt u < tolerance) rnew0

      -- initial conditions for the solver loop
      start :: Acc (Scalar Bool, Scalar Int, Field R, Field R, Field R, Scalar R)
      start = lift (converged0, unit 1, r0, r0, x0, rnew0)

      -- TODO: solver routine (c/linalg.c:223)

      result :: Acc (Scalar Bool, Scalar Int, Field R, Field R, Field R, Scalar R)
      result = error "TODO: solve"
  in
  lift (result ^. _1, result ^. _5)


-- | Compute the inner product of two vectors @x@ and @y@.
--
dot :: Num e => Acc (Vector e) -> Acc (Vector e) -> Acc (Scalar e)
dot xs ys = error "TODO: dot"


-- | Compute the L^2 norm (Euclidean distance) of the given vector
--
norm2 :: Floating e => Acc (Vector e) -> Acc (Scalar e)
norm2 xs = error "TODO: norm2"


-- | The problem is over a rectangular grid of (nx * ny) points. A finite volume
-- discritisation and method of lines gives the following ordinary differential
-- equation for each grid point:
--
-- \[
--    \frac{d s_{ij}}{d t} = \frac{D}{\Delta x^2} [ -4s_{ij} + s_{i-1,j} + s_{i+1,j} + s_{i,j-1} + s_{i,j+1} ] + R s_{ij}(1 - s_{ij})
-- \]
--
-- Which is expressed as the following nonlinear problem:
--
-- \[
--    f_{ij} = [ -(4 + \alpha)s_{ij} + s_{i-1,j} + s_{i+1,j} + s_{i,j-1} + s_{i,j+1} + \beta s_{ij}(1 - s_{ij}) ]^{k+1} + \alpha s_{ij}^k = 0
-- \]
--
-- where s^k is the field at the previous timestep, and s^(k+1) is the current
-- guess of the field at the new time step.
--
-- This corresponds to the 'diffusion' function in 'c/operators.c'. Note that
-- the C implementation inputs the previous grid solution (which I called x0
-- here, and they call X (a macro for x_old)) via a global variable.
--
diffusion
    :: Acc (Scalar R) -- time step
    -> Acc (Field R)  -- solution at time k
    -> Acc (Field R)  -- current approximation of solution at time (k+1)
    -> Acc (Field R)  -- new approximation
diffusion (the -> dt) x0 x1 = error "TODO: diffusion"
  where
    -- Accelerate uses nested tuples to specify stencil patterns, which is
    -- particularly useful for 1D and 2D stencils as they visually represent
    -- which elements of the array are being accessed. Note that we can use as
    -- many or as few of the elements of the stencil as we wish.
    --
    f :: Stencil3x3 R -> Stencil3x3 R -> Exp R
    f s0 s1 = error "TODO: diffusion"

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
    dirichlet :: Boundary (Field R)
    dirichlet = function (\_ -> 0)

    -- simulation constants
    Z :. ny :. nx = unlift (shape x0)
    dx            = 1.0 / (fromIntegral (nx - 1))
    dy            = 1.0 / (fromIntegral (ny - 1))
    dxy           = dx * dy
    alpha         = dxy / (1.0 * dt)            -- diffusion coefficient D = 1.0
    beta          = 1000.0 * dxy

