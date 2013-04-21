module Genetics ( 
                  Configuration (..),
                  defaultConfig,
                  Genetic (mutate, crossover, generate)
                )
                where

import System.Random
import Data.Word
import Data.List
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
   mutate c (BitString bs) = i >>= \n -> return . BitString $ B.append (first n) (flipBit $ second n)
      where
         i              = randomRIO (0,(chromosome_length c)-1)
         first index    = fst $ B.splitAt index bs
         second index   = snd $ B.splitAt index bs
         flipBit :: B.ByteString -> B.ByteString
         flipBit s   | B.head s == 0 = B.cons 1 $ B.tail s
                     | B.head s == 1 = B.cons 0 $ B.tail s
                     | otherwise     = error "FlipBit: Non 0-1 Word8 in ByteString."
         
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
 
runGA :: (Genetic g) => Configuration g -> IO [g]
runGA c = do
               pop <- generate_population c
               breedPop c pop (head pop) 1

generate_population :: (Genetic g) => Configuration g -> IO [g]
generate_population c = helper (chromosomeSize c) []
   where
--    helper :: (Genetic g) => Int -> [g] -> IO [g]
      helper i xs = do
                     if i == 0
                        then
                           return xs
                        else do
                              chromosome <- (generate c)
                              helper (i-1) (chromosome:xs)

breedPop :: (Genetic g) => Configuration g -> [g] -> g -> Int -> IO [g]
breedPop c pop prevMax i = do
                              if (maxGen c) < i
                                 then
                                    do
                                       putStrLn "Maximum generation reached."
                                       return [prevMax]
                                 else
                                    if length filterMatch > 0
                                       then do
                                          putStrLn $ "Target met in "++show i++" generations."
                                          return filterMatch
                                       else do
                                          putStrLn $ "Breeding generation "++show i++"..."
                                          newGen <- newGeneration c tupifiedPop []
                                          breedPop c newGen thisMax $ i+1
   where
      fitnessList = map (fitnessFunction c) pop
      tupifiedPop = zip pop fitnessList
      thisMax     = fst $ maximumBy compareTup (prevMaxTup:tupifiedPop)
      prevMaxTup  = (prevMax, fitnessFunction c $ prevMax)
      filterMatch = [fst x | x <- tupifiedPop, snd x >= (targetFitness c)]
      compareTup t1 t2 = compare (snd t1) (snd t2) 

newGeneration :: (Genetic g) => Configuration g -> [(g,Float)] -> [g] -> IO [g]
newGeneration c old new | populationSize c == length new       = return new
                        | otherwise                            = do
                                                                  (parent1, parent2) <- rouletteSelect c old
                                                                  (child1, child2)   <- checkCrossover c parent1 parent2--crossover
                                                                  (child1', child2') <- checkMutation c child1 child2--mutation
                                                                  if populationSize c == (length new) - 1
                                                                     then
                                                                        newGeneration c old (child1':new)
                                                                     else
                                                                        newGeneration c old (child1':child2':new)

rouletteSelect :: (Genetic g) => Configuration g -> [(g,Float)] -> IO (g,g)
rouletteSelect c xs = child >>= (\w -> child >>= (\y -> return (w,y)))
   where
      rouletteHelper :: (Genetic g) => [(g,Float)] -> Float -> g
      rouletteHelper (x:xs) f | f <= snd x = fst x
                              | otherwise  = rouletteHelper xs (f - snd x)
      child          = (newStdGen >>= return . fst . randomR (0.0, tots_fitness)) >>= return . rouletteHelper xs
      tots_fitness   = foldl (\w -> (\s -> snd s + w)) 0 xs

checkCrossover :: (Genetic g) => Configuration g -> g -> g -> IO (g,g)
checkCrossover c parent1 parent2 = checkFloat (0.0,1.0) (crossoverRate c) >>= \b -> if b then crossover c parent1 parent2 else return (parent1,parent2)

checkMutation ::(Genetic g) => Configuration g -> g -> g -> IO (g,g)
checkMutation c child1 child2 = checkFloat (0.0,1.0) (mutationRate c) >>= \b -> if b then mutate c child1 >>= \a -> mutate c child2 >>= \b -> return (a,b) else return (child1,child2)

checkFloat :: (Float,Float) -> Float -> IO Bool
checkFloat pair f = getRandom pair >>= (\n -> return $ f >= n)

getRandom (lower,upper) = newStdGen >>= return . fst . randomR (lower, upper)








