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
import qualified Data.ByteString as B

type Genetic a where
   mutate      :: Configuration -> a -> IO a
   crossover   :: Configuration -> a -> a -> IO (a,a)
   generate    :: Configuration -> IO a

data BitString = BitString B.ByteString
data OrderedVals = OrderedVals B.ByteString

zero = BitString $ B.singleton 0
one  = BitString $ B.singleton 1

generate_bitstring :: Configuration -> IO BitString
generate_bitstring c = newStdGen >>= (\g -> return . BitString . take (geneSize c) $ randomRs ('0','1') g)

generate_bitstring_chromosome :: Configuration -> IO [BitString]
generate_bitstring_chromosome c = helper (chromosomeSize c) []
   where
      helper :: Int -> [BitString] -> IO [BitString]
      helper i acc = do
                     if (i <= 0)
                        then
                           return acc
                        else do
                           gene <- generate_bitstring c
                           helper (i-1) (gene:acc)

instance Genetic BitString where
   mutate :: Configuration -> BitString -> IO BitString
   mutate c bs = i >>= \n -> return flipBit n bs
      where
         i              = randomRIO (0,(chromosome_length c)-1)
         (head, rest)   = B.splitAt i bs
         flipBit :: BitString -> BitString
         flipBit (Bitstring s)   | B.head s == 0 = BitString $ B.cons 1 $ B.tail s
                                 | B.head s == 1 = BitString $ B.cons 0 $ B.tail s
                                 | otherwise     = error "FlipBit: Non 0-1 Word8 in BitString."
         
   crossover :: Configuration -> BitString -> BitString -> IO (BitString, BitString)
   crossover c parent1 parent2 = return (parent1, parent2)

   generation :: Configuration -> IO BitString
   generation c = word8s >>= \w -> return . BitString $ B.pack w
      where
         word8s :: IO [Word8]
         word8s = newStdGen >>= \g -> return . take . chromosome_length c $ randomRs (0,1) g :: IO [Word8]

chromosome_length :: Configuration -> Int
chromosome_length c = (geneSize c) * (chromosomeSize c)

data Configuration = Config { geneSize             :: Int,
                              chromosomeSize       :: Int,
                              populationSize       :: Int,
                              maxGen               :: Int,
                              targetFitness        :: Float,
                              fitnessFunction      :: Chromosome -> Float,
                              crossoverRate        :: Float,
                              ordered              :: Bool,
                              mutationRate         :: Float
                            }

defaultConfig = Config { geneSize = 0, chromosomeSize = 0, populationSize = 0, maxGen = 0, targetFitness = 0, crossoverRate = 0, mutationRate = 0, ordered = False, fitnessFunction = (\w -> 0.0) }

data Gene = BitString [Char] | Val Int deriving (Show, Eq)
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
                              putStrLn $ "Highest fitness: " ++ show fmaximum
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
                                                                     (child1, child2)   <- checkCrossover c parent1 parent2
                                                                     (child1', child2') <- checkMutation c child1 child2
                                                                     newGeneration c old (child1:new)
                        | otherwise =                             do
                                                                     (parent1, parent2) <- rouletteSelect c old
                                                                     (child1, child2)   <- checkCrossover c parent1 parent2
                                                                     (child1', child2') <- checkMutation c child1 child2
                                                                     newGeneration c old (child1':child2':new)

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

checkCrossover :: Configuration -> Chromosome -> Chromosome -> IO (Chromosome, Chromosome)
checkCrossover c parent1 parent2 = checkFloat (0.0, 1.0) (crossoverRate c) >>= \v -> if v then crossover c parent1 parent2 else return (parent1,parent2)

checkMutation :: Configuration -> Chromosome -> Chromosome -> IO (Chromosome, Chromosome)
checkMutation c child1 child2 = checkFloat (0.0, 1.0) (mutationRate c) >>= \v -> if v then mutate c child1 child2 else (return (child1, child2))

checkFloat :: (Float, Float) -> Float -> IO Bool
checkFloat (lower, upper) p = getRandom (lower, upper) >>= (\n -> return (p >= n))

mutate :: Configuration -> Chromosome -> Chromosome -> IO (Chromosome, Chromosome)
mutate c child1 child2  | not $ ordered c = unorderedMutation c child1 child2
                        | otherwise       = orderedMutation c child1 child2

orderedMutation :: Configuration -> Chromosome -> Chromosome -> IO (Chromosome, Chromosome)
orderedMutation c child1 child2 = swapValues >>= \l -> return $ (swap (head l) (head . tail $ l) child1, swap (head . tail . tail $ l) (last l) child2)
   where
      swapValues = sequence [getRandom (0, cSize), getRandom (0, cSize), getRandom (0, cSize), getRandom (0, cSize)]
      cSize = (chromosomeSize c) - 1 
      swap :: Int -> Int -> Chromosome -> Chromosome
      swap x y child = part1 ++ (val_y:part2) ++ (val_x:part3)
         where
            x' = min x y
            y' = max x y
            val_x = child !! x'
            val_y = child !! y'
            part1 = take x' child
            part2 = take (y' - x' - 1) $ drop (x'+1) child 
            part3 = drop (y'+1) child

unorderedMutation :: Configuration -> Chromosome -> Chromosome -> IO (Chromosome, Chromosome)
unorderedMutation c child1 child2 = getRandom (1, numberOfBitsInChromosome) >>= (\x -> getRandom (0, numberOfBitsInChromosome) >>= (\y -> return (reformChromosome (geneSize c) . mutateHelper x $ collapseChromosome (geneSize c) child1, reformChromosome (geneSize c) . mutateHelper y $ collapseChromosome (geneSize c) child2)))
   where
      numberOfBitsInChromosome = (geneSize c) * (chromosomeSize c)
      flipBit :: Char -> Char
      flipBit c = if c == '1' then '0' else '1'
      mutateHelper i []     = []
      mutateHelper i (x:xs) = if i == 1 then (flipBit x):(mutateHelper (i-1) xs) else x:(mutateHelper (i-1) xs)

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

generate_bitstring :: Configuration -> IO BitString
generate_bitstring c = newStdGen >>= (\g -> return . BitString . take (geneSize c) $ randomRs ('0','1') g)

generate_bitstring_chromosome :: Configuration -> IO [BitString]
generate_bitstring_chromosome c = helper (chromosomeSize c) []
   where
      helper :: Int -> [BitString] -> IO [BitString]
      helper i acc = do
                     if (i <= 0)
                        then
                           return acc
                        else do
                           gene <- generate_bitstring c
                           helper (i-1) (gene:acc)


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


