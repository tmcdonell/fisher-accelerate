{-# LANGUAGE ConstraintKinds  #-}
{-# LANGUAGE FlexibleContexts #-}

-- A small program which solves the 2D Fisher equation using second-order finite
-- differences.
--

module Main where

-- Accelerate supports multiple backends. Change the input here to switch which
-- backend to use.
--
-- import Data.Array.Accelerate.LLVM.Native                            as CPU
import Data.Array.Accelerate.LLVM.PTX                               as PTX

-- For the exercise, you will need to implement functions in this module.
--
import Solve


-- Other imports for this module
import Data.Array.Accelerate                                        as A
import Data.Array.Accelerate.Array.Data                             ( ptrsOfArrayData )
import Data.Array.Accelerate.Array.Sugar                            ( Array(..) )
import Data.Array.Accelerate.Debug                                  ( showFFloatSIBase )
import qualified Data.Array.Accelerate.Debug                        as Debug

import Control.Exception
import Control.Monad
import Data.Time.Clock
import Foreign.Storable
import System.CPUTime
import System.Environment
import System.IO
import Text.Printf
import Prelude                                                      as P


-- The main simulation loop
--
main :: IO ()
main = do
  -- This will strip out the "+ACC ... -ACC" from the command line arguments,
  -- which makes handling the input parameters a bit simpler. Alternatively,
  -- they can be specified using the ACCELERATE_FLAGS environment variable.
  --
  -- You may want to turn on a few of the debugging options to get some more
  -- information about what the compiler is doing.
  --
  Debug.accInit

  -- The command line arguments, stripped of Accelerate RTS options
  argv <- getArgs

  let
      dx                    = 1 / (P.fromIntegral nx - 1) :: R
      (nx, ny, steps, tend) =
        case argv of
          [a,b,c,d]
            | [(a',[])] <- reads a
            , [(b',[])] <- reads b
            , [(c',[])] <- reads c
            , [(d',[])] <- reads d -> (a',b',c',d')
          _                        -> usage

      -- The simulation timestep
      --
      -- This is put into a Scalar (single element) array so that our main loop
      -- can be put into the form of a function from arrays to arrays; i.e.
      -- (Acc a -> Acc b -> ... Acc c).
      --
      dt :: Scalar R
      dt = fromList Z [tend / P.fromIntegral steps]

      -- Generate the initial conditions
      --
      x0 :: Field R
      x0 = run $ initialise nx ny

      -- Compile and run the main loop
      --
      -- There are different ways that this could be done. Experiment a bit
      -- here...
      --
      go i x
        | i P.>= steps  = x
        | otherwise     =
            let (ok, x') = run $ step (use dt) (use x)
            in if indexArray ok Z
                 then go (i+1) x'
                 else error (printf "step %d error: non-linear iterations failed to converge" i)

  -- Start EKG monitor
  --
  Debug.beginMonitoring

  -- Show the welcome banner
  --
  nx `seq` ny `seq` steps `seq` tend `seq` return ()  -- if there are errors, get them early
  printf "========================================================================\n"
  printf "                      Welcome to mini-stencil!\n"
  printf "mesh :: %d * %d, dx = %.4f\n" nx ny dx
  printf "time :: %d, time steps from 0 .. %f\n" steps tend
  printf "========================================================================\n"

  -- Run the simulation and store the final timestep result to file. You can
  -- open the 'output.bov' file in a program like VisIt to check that it looks
  -- right.
  --
  r   <- elapsed (evaluate (go 0 x0))
  writeBOV steps tend r



-- Set the initial conditions: a circle of concentration of 0.1 centred at
-- (1/4, 1/4) of a unit grid, with radius 1/8 units.
--
initialise
    :: Int
    -> Int
    -> Acc (Field R)
initialise xdim ydim =
  let
      dx  = 1 / (P.fromIntegral xdim - 1 :: R)
      dy  = 1 / (P.fromIntegral ydim - 1 :: R)
      r   = 1 / 8

      f :: Exp DIM2 -> Exp R
      f ix  =
        let Z :. j :. i = unlift ix
            x           = A.fromIntegral (i - 1) * A.constant dx
            y           = A.fromIntegral (j - 1) * A.constant dy
        in
        (x - 0.25) * (x - 0.25) + (y - 0.25) * (y - 0.25) A.< A.constant (r*r) ? (0.1, 0.0)
  in
  generate (constant (Z:.ydim:.xdim)) f


-- Save the result of the given simulation step to file for later visualisation.
-- The file can be opened with the VisIt program. The script 'tools/phi.py' can
-- be used to generate the picture (via VisIt) as well.
--
-- https://wci.llnl.gov/simulation/computer-codes/visit
--
writeBOV :: Int -> R -> Field R -> IO ()
writeBOV _ t arr@(Array _ adata) = do
  let
      outmeta       = "output.bov"
      outdata       = "output.bin"
      --
      Z :. ny :. nx = arrayShape arr

  -- add file metadata
  withFile outmeta AppendMode $ \h -> do
    hPutStrLn h (printf "TIME: %f"              t)
    hPutStrLn h (printf "DATA_FILE: %s"         outdata)
    hPutStrLn h (printf "DATA_SIZE: %d, %d, 1"  nx ny)
    hPutStrLn h (printf "DATA_FORMAT: %s"       (if sizeOf (undefined::R) P.== 4 then "FLOAT" else "DOUBLE"))
    hPutStrLn h (printf "VARIABLE: phi");
    hPutStrLn h (printf "CENTERING: nodal")
    hPutStrLn h (printf "BRICK_SIZE: 1.0 1.0 1.0")
    -- hPutStrLn h (printf "DATA_ENDIAN: LITTLE")
    -- hPutStrLn h (printf "BYTE_OFFSET: 4")

  -- dump the data
  withFile outdata WriteMode $ \h ->
    hPutBuf h (ptrsOfArrayData adata) (arraySize (arrayShape arr) * sizeOf (undefined::R))


usage :: a
usage = error $ unlines
  [ "fisher-accelerate"
  , ""
  , "Usage: fisher-accelerate nx ny s t"
  , "  nx   number of grid points in the x-direction"
  , "  ny   number of grid points in the y-direction"
  , "  s    number of simulation steps"
  , "  t    total simulation time"
  , ""
  , "Example parameters:"
  , "  128 128 100 0.0025"
  , "  256 256 100 0.01"
  ]

elapsed :: IO a -> IO a
elapsed action = do
  wall0 <- getCurrentTime
  cpu0  <- getCPUTime
  res   <- action
  wall1 <- getCurrentTime
  cpu1  <- getCPUTime
  --
  let wallTime = realToFrac (diffUTCTime wall1 wall0)
      cpuTime  = P.fromIntegral (cpu1 - cpu0) * 1E-12
  --
  printf "Simulation took %s (wall), %s (cpu), %.2f x speedup\n"
    (showFFloatSIBase (Just 3) 1000 wallTime "s")
    (showFFloatSIBase (Just 3) 1000 cpuTime  "s")
    (cpuTime / wallTime :: Double)
  --
  return res

