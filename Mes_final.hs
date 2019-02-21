-- To włącza pewne rozszerzenie składni, które czyni korzystanie z Haskellowych
-- rekordów sporo wygodniejszym w pewnych sytuacjach.
-- https://ocharles.org.uk/posts/2014-12-04-record-wildcards.html
{-# LANGUAGE RecordWildCards #-}
module Main where

import Numeric.LinearAlgebra
import Graphics.EasyPlot as GEP
import Control.Applicative
import Control.Monad

type Point = Double
type Val = Double
type Fun = Point -> Val
type DoF = Int

data Interval = Interval {
    leftPoint  :: Point,
    rightPoint :: Point
}

intervalLength :: Interval -> Val
intervalLength (Interval a b) = b - a

inside :: Interval -> Point -> Bool
inside (Interval a b) x = a <= x && x <= b

data Mesh = Mesh {
    domain   :: Interval,
    elements :: Int
}

elementSize :: Mesh -> Val
elementSize (Mesh domain n) = intervalLength domain / fromIntegral n

meshPoint :: Mesh -> Int -> Point
meshPoint mesh k = a + h * fromIntegral k
    where a = leftPoint $ domain mesh
          h = elementSize mesh

element :: Mesh -> Int -> Interval
element mesh k = Interval left right
    where left  = meshPoint mesh k
          right = meshPoint mesh (k + 1)


basisFunction :: Mesh -> DoF -> Fun
basisFunction mesh i x
    | inElem (i - 1) = s * (x - x0)
    | inElem i       = s * (x1 - x)
    | otherwise      = 0
    where inElem i = inside (element mesh i) x
          s  = 1 / elementSize mesh
          x0 = meshPoint mesh (i - 1)
          x1 = meshPoint mesh (i + 1)

basisFunctionDer :: Mesh -> DoF -> Fun
basisFunctionDer mesh i x
    | inElem (i - 1) = s
    | inElem i       = -s
    | otherwise      = 0
    where inElem i = inside (element mesh i) x
          s = 1 / elementSize mesh

integrate' :: Interval -> Int -> Fun -> Val
integrate' iv@(Interval a b) n f = sum [ w * (f (x k)) | k <- [0..n] ]
    where w   = intervalLength iv / fromIntegral n
          x k = a + w * fromIntegral k

data Problem = Problem {
    a, b, c, f :: Fun, -- równanie
    u0 :: Val,            -- warunek Dirichleta
    beta, gamma :: Val   -- warunek Robina
}



infixl 6 .+, .-
infixl 7 .*


(.+), (.*), (.-) :: Fun -> Fun -> Fun
(.+) = liftA2 (+)
(.*) = liftA2 (*)
(.-) = liftA2 (-)

quadPoints :: Int
quadPoints = 10000

formB :: Problem -> Mesh -> Fun -> Fun -> Fun -> Fun -> Val
formB Problem{..} mesh u u' v v' =
    int (b.*u'.*v .- a.*u'.*v' .+ c.*u.*v) - beta * (u 0) * (v 0)
    where int = integrate' (domain mesh) quadPoints

formL :: Problem -> Mesh -> Fun -> Val
formL Problem{..} mesh v =
    int (f .* v) - gamma * (v 0)
    where int = integrate' (domain mesh) quadPoints

linearCombination :: Mesh -> Vector Double -> Fun
linearCombination mesh coeffs x = dot coeffs (build n val)
    where val i = basisFunction mesh (round i) x
          n = elements mesh

solve :: Problem -> Mesh -> Fun
solve prob@Problem{..} mesh = shift .+ linearCombination mesh coeffs
    where coeffs = matB <\> vecL
          matB = build (n, n) entryB :: Matrix Double
          vecL = build n      entryL :: Vector Double
          entryB i j = b (u j) (u' j) (u i) (u' i)
          entryL i   = l (u i) - b shift shift' (u i) (u' i)
          b = formB prob mesh
          l = formL prob mesh
          shift  = const u0
          shift' = const 0
          u  i = basisFunction    mesh (round i)
          u' i = basisFunctionDer mesh (round i)
          n = elements mesh

plotFunction :: Fun -> Interval -> IO ()
plotFunction f (Interval a b) = void $ plot X11 $ Function2D [Title "Solution"] [GEP.Range a b] f

mesh = Mesh {
    domain = Interval 0 1,
    elements = 10
}

problem = Problem {
    a = \x -> x^2 + 1,
    b = \x -> - 4 * x,
    c = \x -> x + 1,
    f = \x -> -x^3 + 2 * x^2 - x - 2,

    u0 = 0,

    gamma = -1, beta = 0
}

doStuff :: IO ()
doStuff = plotFunction (solve problem mesh) (domain mesh)

main = do
 doStuff
 return 0
