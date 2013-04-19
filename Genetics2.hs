module Genetics ( 
                  Configuration (..),
                  defaultConfig,
                  Genetic (mutate, crossover, generate)
                )
                where

import System.Random
import Data.Word
import qualified Data.ByteString as B

class Genetic a where
   mutate      :: Configuration a -> a -> IO a
   crossover   :: Configuration a -> a -> a -> IO (a,a)
   generate    :: Configuration a -> IO a

data BitString = BitString B.ByteString 
data OrderedVals = OrderedVals B.ByteString 

--zero = BitString $ B.singleton 0
--one  = BitString $ B.singleton 1

instance Genetic BitString where

   --mutate :: Configuration BitString -> BitString -> IO BitString
   mutate c (BitString bs) = i >>= \n -> return . BitString $ B.append (first n) (second n)
      where
         i              = randomRIO (0,(chromosome_length c)-1)
         first index    = fst $ B.splitAt index bs
         second index   = snd $ B.splitAt index bs
         flipBit :: BitString -> BitString
         flipBit (BitString s)   | B.head s == 0 = BitString $ B.cons 1 $ B.tail s
                                 | B.head s == 1 = BitString $ B.cons 0 $ B.tail s
                                 | otherwise     = error "FlipBit: Non 0-1 Word8 in BitString."
         
   --crossover :: Configuration BitString -> BitString -> BitString -> IO (BitString, BitString)
   crossover c parent1 parent2 = i >>= \n -> return $ onePointCrossover n parent1 parent2
      where
         i = randomRIO (0, (chromosome_length c)-1)
         onePointCrossover :: Int -> BitString -> BitString -> (BitString, BitString)
         onePointCrossover i (BitString p1) (BitString p2)  = let ((head1, tail2),(head2,tail1)) = ((B.splitAt i p1),(B.splitAt i p2)) in ((BitString $ B.append head1 tail1), (BitString $ B.append head2 tail2))

   --generate :: Configuration BitString -> IO BitString
   generate c = bits >>= \w -> return . BitString $ B.pack w
      where
         bits :: IO [Word8]
         bits = newStdGen >>= \g -> return . take (chromosome_length c) $ randomRs (0,1) g :: IO [Word8]

instance Show BitString where
--show :: BitString -> String
   show (BitString bs) = B.foldr (helper) "" bs
      where
         helper word acc   | word == 0 = '0':acc
                           | word == 1 = '1':acc
                           | otherwise = error "Show: Non 0-1 bytes in BitString."

chromosome_length :: Configuration a -> Int
chromosome_length c = (geneSize c) * (chromosomeSize c)

data Configuration a = Config { 
                              geneSize             :: Int,
                              chromosomeSize       :: Int,
                              populationSize       :: Int,
                              maxGen               :: Int,
                              targetFitness        :: Float,
                              fitnessFunction      :: a -> Float,
                              crossoverRate        :: Float,
                              ordered              :: Bool,
                              mutationRate         :: Float
                              }

defaultConfig = Config { geneSize = 0, chromosomeSize = 0, populationSize = 0, maxGen = 0, targetFitness = 0, crossoverRate = 0, mutationRate = 0, ordered = False, fitnessFunction = (\w -> 0.0) }

--data Gene = BitString [Char] | Val Int deriving (Show, Eq)
--type Chromosome = [Gene]
--type Population = [Chromosome]
 
--runGA :: Configuration -> IO ()
--runGA c = do
--               pop <- generate_population c
--               breedPop c pop 1

--generate_population :: (Genetic g) => Configuration -> IO [g]
--generate_population
