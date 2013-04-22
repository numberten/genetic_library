import Genetics2
import Text.Regex
import qualified Data.ByteString as B

config = defaultConfig {   populationSize = 100,
                           fitnessFunction = f,
                           geneSize = 4,
                           chromosomeSize = 9,
                           targetFitness = 50.0,
                           maxGen = 100,
                           crossoverRate = 0.8,
                           mutationRate = 0.1,
                           preciseFitness = True
                       }

f :: BitString -> Float
f bs = if fit == 50 then 50 else evaluateFitness fit
   where
      fit = evaluateChromosome bs


evaluateFitness :: Float -> Float
evaluateFitness v
		| or [show v == "Infinity", show v == "NaN"] = 0 --Invalid chromosomes. Infinity/NaN, default 0 fitness.
		| otherwise = abs $ 1 / ((targetFitness config) - v)	 	

-- ###################################### CHROMOSOME DECODING CODE BELOW ########################################### 
decodeChromosome :: BitString -> [Char]
decodeChromosome (BitString xs) 
		| mod (B.length xs) 4 /= 0 = error "Chromosome has corrupted genes."
		| (show . BitString $ B.take 4 xs) == "0000" = '0':(decodeChromosome (BitString $ B.drop 4 xs))
		| (show . BitString $ B.take 4 xs) == "0001" = '1':(decodeChromosome (BitString $ B.drop 4 xs))
		| (show . BitString $ B.take 4 xs) == "0010" = '2':(decodeChromosome (BitString $ B.drop 4 xs))
		| (show . BitString $ B.take 4 xs) == "0011" = '3':(decodeChromosome (BitString $ B.drop 4 xs))
		| (show . BitString $ B.take 4 xs) == "0100" = '4':(decodeChromosome (BitString $ B.drop 4 xs))
		| (show . BitString $ B.take 4 xs) == "0101" = '5':(decodeChromosome (BitString $ B.drop 4 xs))
		| (show . BitString $ B.take 4 xs) == "0110" = '6':(decodeChromosome (BitString $ B.drop 4 xs))
		| (show . BitString $ B.take 4 xs) == "0111" = '7':(decodeChromosome (BitString $ B.drop 4 xs))
		| (show . BitString $ B.take 4 xs) == "1000" = '8':(decodeChromosome (BitString $ B.drop 4 xs))
		| (show . BitString $ B.take 4 xs) == "1001" = '9':(decodeChromosome (BitString $ B.drop 4 xs))
		| (show . BitString $ B.take 4 xs) == "1010" = '+':(decodeChromosome (BitString $ B.drop 4 xs))
		| (show . BitString $ B.take 4 xs) == "1011" = '-':(decodeChromosome (BitString $ B.drop 4 xs))
		| (show . BitString $ B.take 4 xs) == "1100" = '*':(decodeChromosome (BitString $ B.drop 4 xs))
		| (show . BitString $ B.take 4 xs) == "1101" = '/':(decodeChromosome (BitString $ B.drop 4 xs))
		| (show . BitString $ B.take 4 xs) == "1110" = ' ':(decodeChromosome (BitString $ B.drop 4 xs))
		| (show . BitString $ B.take 4 xs) == "1111" = ' ':(decodeChromosome (BitString $ B.drop 4 xs))
                | otherwise = []

spaceOperators :: [Char] -> [Char]
spaceOperators = subMinus . subPlus . subTimes. subDivide 

subMinus :: [Char] -> [Char]
subMinus xs = subRegex (mkRegex "\\-") xs " - "

subPlus :: [Char] -> [Char]
subPlus xs = subRegex (mkRegex "\\+") xs " + "

subTimes :: [Char] -> [Char]
subTimes xs = subRegex (mkRegex "\\*") xs " * "

subDivide :: [Char] -> [Char]
subDivide xs = subRegex (mkRegex "\\/") xs " / "

subSpace :: [Char] -> [Char]
subSpace xs = subRegex (mkRegex " ") xs ""

operator :: Char -> Bool
operator xs = or [xs == '+', xs == '-', xs == '/', xs == '*']

removeOpHeads :: [Char] -> [Char]
removeOpHeads [] = []
removeOpHeads z@(x:xs)
		| operator x = removeOpHeads xs
		| otherwise  = z

removeOpLasts :: [Char] -> [Char]
removeOpLasts [] = []
removeOpLasts xs 
		| operator $ last xs = removeOpLasts $ init xs
		| otherwise = xs

removeOpDoubles :: [Char] -> [Char]
removeOpDoubles []  = []
removeOpDoubles [x] = [x]
removeOpDoubles (x:y:ys)
		| and [operator x, operator y] = (removeOpDoubles $ x:ys)
		| otherwise = x:(removeOpDoubles $ y:ys)

--read . head . [[Char]]
evaluateExpression :: [Char] -> Float
evaluateExpression = read . head . reverse . (:) "Infinity" . reverse . foldl helperFunction [] . words . spaceOperators . removeOpDoubles . removeOpHeads . removeOpLasts . subSpace
	where	helperFunction ("+":y:ys) n = (show $ (read n :: Float) + (read y :: Float)):ys
		helperFunction ("-":y:ys) n = (show $ (read y :: Float) - (read n :: Float)):ys
		helperFunction ("*":y:ys) n = (show $ (read n :: Float) * (read y :: Float)):ys
		helperFunction ("/":y:ys) n = (show $ (read y :: Float) / (read n :: Float)):ys
		helperFunction xs n = n:xs

evaluateChromosome :: BitString -> Float
evaluateChromosome = evaluateExpression . decodeChromosome 

-- ###################################### CHROMOSOME DECODING CODE ABOVE ########################################### 






decodeChromosomeC [] = []
decodeChromosomeC xs
		| mod (length xs) 4 /= 0 = error "Chromosome has corrupted genes."
		| take 4 xs == "0000" = '0':(decodeChromosomeC $ drop 4 xs)
		| take 4 xs == "0001" = '1':(decodeChromosomeC $ drop 4 xs)
		| take 4 xs == "0010" = '2':(decodeChromosomeC $ drop 4 xs)
		| take 4 xs == "0011" = '3':(decodeChromosomeC $ drop 4 xs)
		| take 4 xs == "0100" = '4':(decodeChromosomeC $ drop 4 xs)
		| take 4 xs == "0101" = '5':(decodeChromosomeC $ drop 4 xs)
		| take 4 xs == "0110" = '6':(decodeChromosomeC $ drop 4 xs)
		| take 4 xs == "0111" = '7':(decodeChromosomeC $ drop 4 xs)
		| take 4 xs == "1000" = '8':(decodeChromosomeC $ drop 4 xs)
		| take 4 xs == "1001" = '9':(decodeChromosomeC $ drop 4 xs)
		| take 4 xs == "1010" = '+':(decodeChromosomeC $ drop 4 xs)
		| take 4 xs == "1011" = '-':(decodeChromosomeC $ drop 4 xs)
		| take 4 xs == "1100" = '*':(decodeChromosomeC $ drop 4 xs)
		| take 4 xs == "1101" = '/':(decodeChromosomeC $ drop 4 xs)
		| take 4 xs == "1110" = ' ':(decodeChromosomeC $ drop 4 xs)
		| take 4 xs == "1111" = ' ':(decodeChromosomeC $ drop 4 xs)
