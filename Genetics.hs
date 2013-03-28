module Genetics ( 
                  Configuration (..),
                  runGA
                )
                where

import System.Random

data Configuration = Config { geneSize             :: Int,
                              chromosomeSize       :: Int,
                              populationSize       :: Int,
                              maxGen               :: Int,
                              targetFitness        :: Float,
                              fitnessFunction      :: Chromosome -> Float,
                              crossoverRate        :: Float
                            }

type Gene = [Char]
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
                        putStrLn ("Target met in generation " ++ show i ++ ".\n" ++ (showTargets targetsc))
                     else
                        if (maxGen c < i)
                           then
                              putStrLn "Maximum generation reached."
                           else do
                              putStrLn ("Finished with iteration " ++ show i ++ ".\nBreeding population...")
                              newGen <- newGeneration c p []
                              breedPop c newGen (i+1)
                  where 
                  fitnessList = map (fitnessFunction c) p
                  targetsc = [x | x <- p, (fitnessFunction c) x == targetFitness c]
                  targets = filter (== targetFitness c) fitnessList
                  fmaximum = maximum fitnessList
                  faverage = ftotal / (fromIntegral $ length p)
                  ftotal = sum fitnessList

newGeneration :: Configuration -> Population -> Population -> IO Population
newGeneration c old new | populationSize c == length new = return new
                        | populationSize c == (length new) - 1 =  do
                                                                     (parent1, parent2) <- rouletteSelect c old
                                                                     (child1, _)        <- crossover c parent1 parent2
                                                                     newGeneration c old (child1:new)
                        | otherwise =                             do
                                                                     (parent1, parent2) <- rouletteSelect c old
                                                                     (child1, child2)   <- crossover c parent1 parent2
                                                                     newGeneration c old (child1:child2:new)

showTargets :: [Chromosome] -> String
showTargets []     = ""
showTargets [x]    = "Chromosome: "++ (concat x) ++ "."
showTargets (x:xs) = "Chromosome: "++ (concat x) ++ ","

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

-- crossover :
crossover :: Configuration -> Chromosome -> Chromosome -> IO (Chromosome,Chromosome)
crossover c c1 c2 = (newStdGen >>= return . fst . randomR (0.0, 1.0)) >>= (\x -> if (x > crossoverRate c) then return (c1, c2) else getFlipBit c >>= (\z -> onePointCrossover (geneSize c) z c1 c2))

--returns index of random bit in a Chromosome
getFlipBit :: Configuration -> IO Int
getFlipBit c = newStdGen >>= return . fst . randomR (1, (geneSize c) * (chromosomeSize c))

onePointCrossover :: Int -> Int -> Chromosome -> Chromosome -> IO (Chromosome, Chromosome)
onePointCrossover genesize i c1 c2 = return (breakByGene genesize (take i c1' ++ drop i c2') [], breakByGene genesize (take i c2' ++ drop i c1') [])
               where
                  c1' = concat c1
                  c2' = concat c2
                  breakByGene :: Int -> Gene -> Chromosome -> Chromosome
                  breakByGene _ [] acc = acc
                  breakByGene i cs acc = breakByGene i (drop i cs) ((take i cs):acc)

-- Generate Functions --

generate_gene :: Configuration -> IO Gene
generate_gene c = newStdGen >>= (\g -> return . take (geneSize c) $ randomRs ('0','1') g)

generate_chromosome :: Configuration -> IO Chromosome
generate_chromosome c = w2 (chromosomeSize c) []
   where
      w2 :: Int -> Chromosome -> IO Chromosome
      w2 i acc = do
                     if (i <= 0)
                        then
                           return acc
                        else do
                           gene <- generate_gene c
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
                           w3 (i-1) (chromosome:acc)

