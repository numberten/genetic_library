module Genetics ( 
                  Configuration (..),
                  defaultConfig,
                  Gene(Val),
                  Chromosome,
                  runGA
                )
                where

import System.Random
import Data.List

data Configuration = Config { geneSize             :: Int,
                              chromosomeSize       :: Int,
                              populationSize       :: Int,
                              maxGen               :: Int,
                              targetFitness        :: Float,
                              fitnessFunction      :: Chromosome -> Float,
                              crossoverRate        :: Float,
                              ordered              :: Bool
                            }

defaultConfig = Config { geneSize = 0, chromosomeSize = 0, populationSize = 0, maxGen = 0, targetFitness = 0, crossoverRate = 0, ordered = False, fitnessFunction = (\w -> 0.0) }

data Gene = BitString [Char] | Val Int deriving (Show, Eq)
--type Gene = [Char]
type Chromosome = [Gene]
type Population = [Chromosome]
 
runGA :: Configuration -> IO ()
runGA c = do
               pop <- generate_population c
               breedPop c pop 1

breedPop :: Configuration -> Population -> Int -> IO ()
breedPop c p i = do
                  if (targets /= [])
                     then
                        putStrLn ("Target met in generation " ++ show i ++ ".\n" ++ (showTargets (ordered c) targetsc))
                     else
                        if (maxGen c < i)
                           then
                              putStrLn ("Maximum generation reached.\nMost fit chromosome: " ++ (showTargets (ordered c) mostFitc) ++ " with fitness " ++ show fmaximum)
                           else do
                              putStrLn ("Finished with iteration " ++ show i ++ ".\nBreeding population...")
                              newGen <- newGeneration c p []
                              putStrLn $ "Total fitness: " ++ show ftotal
                              breedPop c newGen (i+1)
                  where 
                  fitnessList = map (fitnessFunction c) p
                  targetsc = [x | x <- p, (fitnessFunction c) x >= targetFitness c]
                  mostFitc = [x | x <- p, (fitnessFunction c) x == fmaximum]
                  targets = filter (>= targetFitness c) fitnessList
                  mostFit = filter (== fmaximum) fitnessList
                  fmaximum = maximum fitnessList
                  faverage = ftotal / (fromIntegral $ length p)
                  ftotal = sum fitnessList

newGeneration :: Configuration -> Population -> Population -> IO Population
newGeneration c old new | populationSize c == length new = return new
                        | populationSize c == (length new) - 1 =  do
                                                                     (parent1, parent2) <- rouletteSelect c old
                                                                     bool <- checkCrossover c
                                                                     if (bool)
                                                                        then do
                                                                           (child1, _) <- crossover c parent1 parent2
                                                                           newGeneration c old (child1:new)
                                                                        else
                                                                           newGeneration c old (parent1:new)
                        | otherwise =                             do
                                                                     (parent1, parent2) <- rouletteSelect c old
                                                                     bool <- checkCrossover c
                                                                     if (bool)
                                                                        then do
                                                                           (child1, child2) <- crossover c parent1 parent2
                                                                           newGeneration c old (child1:child2:new)
                                                                        else
                                                                           newGeneration c old (parent1:parent2:new)

showTargets :: Bool -> [Chromosome] -> String
showTargets ord []     = ""
showTargets ord [x]    | not ord      = "Chromosome: "++ (showUnorderedChromosome x)
                       | otherwise = "Chromosome: "++ (showOrderedChromosome x)
showTargets ord (x:xs) | not ord      = "Chromosome: "++ (showUnorderedChromosome x) ++ "\n" ++ (showTargets ord xs)
                       | otherwise = "Chromosome: "++ (showOrderedChromosome x) ++ "\n" ++ (showTargets ord xs)
   --where
showUnorderedChromosome :: Chromosome -> String
showUnorderedChromosome []                 = ""
showUnorderedChromosome [(BitString x)]    = x ++ "."
showUnorderedChromosome ((BitString x):xs) = x ++ (showUnorderedChromosome xs)
   --where
showOrderedChromosome :: Chromosome -> String
showOrderedChromosome []           = ""
showOrderedChromosome [(Val x)]    = show x ++ "."
showOrderedChromosome ((Val x):xs) = show x ++ ", " ++ (showOrderedChromosome xs)

rouletteSelect :: Configuration -> Population -> IO (Chromosome, Chromosome)
rouletteSelect c p = child1 >>= (\w -> child2 >>= (\y -> return (w,y)))
            where
               child1       = (newStdGen >>= return . fst . randomR (0.0, tots_fitness)) >>= return . rouletteHelper z
               child2       = (newStdGen >>= return . fst . randomR (0.0, tots_fitness)) >>= return . rouletteHelper z
               z            = tupify p (fitnessFunction c)
               tots_fitness = foldl (\w -> (\s -> snd s + w)) 0 z

rouletteHelper :: [(Chromosome, Float)] -> Float -> Chromosome
rouletteHelper (x:xs) f | f <= snd x = fst x
                        | otherwise  = rouletteHelper xs (f - snd x)

tupify :: Population -> (Chromosome -> Float) -> [(Chromosome, Float)]
tupify p f = map (\w -> (w,f w)) p

checkCrossover :: Configuration -> IO Bool
checkCrossover c = getRandom (0.0, 1.0) >>= (\n -> return (crossoverRate c >= n))

crossover :: Configuration -> Chromosome -> Chromosome -> IO (Chromosome,Chromosome)
crossover c parent1 parent2   | not $ ordered c = onePointCrossover c parent1 parent2
                              | otherwise       = crossoverPoint1 >>= (\cp1 -> crossoverPoint2 >>= (\cp2 -> orderedCrossover cp1 cp2 parent1 parent2))
   where
      crossoverPoint1 = getRandom (0, chromosomeSize c)
      crossoverPoint2 = getRandom (0, chromosomeSize c)

onePointCrossover :: Configuration -> Chromosome -> Chromosome -> IO (Chromosome, Chromosome)
onePointCrossover c parent1 parent2 = crossoverPoint >>= (\p -> return (reformChromosome (geneSize c) (take p collapsedParent1 ++ drop p collapsedParent2), reformChromosome (geneSize c) (take p collapsedParent2 ++ drop p collapsedParent1)))
   where
      crossoverPoint   = getRandom (1, (geneSize c) * (chromosomeSize c))
      collapsedParent1 = collapseChromosome (geneSize c) parent1
      collapsedParent2 = collapseChromosome (geneSize c) parent2
      collapseChromosome :: Int -> Chromosome -> [Char]
      collapseChromosome _ []                        = []
      collapseChromosome genesize ((BitString x):xs) = x ++ (collapseChromosome genesize xs)
      reformChromosome :: Int -> [Char] -> Chromosome 
      reformChromosome _ []        = []
      reformChromosome genesize xs = (BitString (take genesize xs)):(reformChromosome genesize (drop genesize xs))

orderedCrossover :: Int -> Int -> Chromosome -> Chromosome -> IO (Chromosome,Chromosome)
orderedCrossover crossoverPoint1 crossoverPoint2 parent1 parent2 = return (child1, child2)
   where
      lowestPoint     = min crossoverPoint1 crossoverPoint2
      highestPoint    = max crossoverPoint1 crossoverPoint2
      firstThirdP1    = take lowestPoint parent1
      secondThirdP1   = take (highestPoint - lowestPoint) (drop lowestPoint parent1)
      thirdThirdP1    = drop highestPoint parent1
      firstThirdP2    = take lowestPoint parent2
      secondThirdP2   = take (highestPoint - lowestPoint) (drop lowestPoint parent2)
      thirdThirdP2    = drop highestPoint parent2
      (l1, l2, l3)    = (length firstThirdP1, length secondThirdP1, length thirdThirdP1)
      child1          = (take l1 remainderC1) ++ secondThirdP1 ++ (drop l1 remainderC1)
      child2          = (take l1 remainderC2) ++ secondThirdP2 ++ (drop l1 remainderC2)
      remainderC1     = (thirdThirdP2 ++ firstThirdP2 ++ secondThirdP2) \\ secondThirdP1
      remainderC2     = (thirdThirdP1 ++ firstThirdP1 ++ secondThirdP1) \\ secondThirdP2

getRandom (lower,upper) = newStdGen >>= return . fst . randomR (lower, upper)

-- Generate Functions --

generate_gene :: Int -> Configuration -> IO Gene
generate_gene i c | not $ ordered c = newStdGen >>= (\g -> return . BitString . take (geneSize c) $ randomRs ('0','1') g)
                  | otherwise       = return $ Val i

generate_chromosome :: Configuration -> IO Chromosome
generate_chromosome c = w2 (chromosomeSize c) []
   where
      w2 :: Int -> Chromosome -> IO Chromosome
      w2 i acc = do
                     if (i <= 0)
                        then
                           return acc
                        else do
                           gene <- generate_gene i c
                           w2 (i-1) (gene:acc)

generate_population :: Configuration -> IO Population
generate_population c = w3 (populationSize c) []
   where
      w3 :: Int -> Population -> IO Population
      w3 i acc = do
                     if (i <= 0)
                        then
                           return acc
                        else do
                           chromosome <- generate_chromosome c
                           shuffled_chromosome <- shuffle chromosome []
                           w3 (i-1) (shuffled_chromosome:acc)

shuffle :: (Eq a) => [a] -> [a] -> IO [a]
shuffle old new = if (old == []) then return new else do
                                                         r <- getRandom (0, length old - 1)
                                                         shuffle (delete (old !! r) old) ((old !! r):new)


