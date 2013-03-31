import Genetics

config = defaultConfig {   chromosomeSize = 10,
                           populationSize = 200,
                           maxGen = 20,
                           targetFitness = 15000,
                           crossoverRate = 0.8,
                           ordered = True,
                           fitnessFunction = f }

f :: Chromosome -> Float
--f = totalDistance
f c = (/) 1.0 . totalDistance (head c) $ c

totalDistance :: Gene -> Chromosome -> Float
totalDistance _ []                          = 0
totalDistance (Val o) [Val i]               = fromIntegral (distanceLookup i o distances)
totalDistance origin z@((Val i):(Val j):xs) = fromIntegral (distanceLookup i j distances) + (totalDistance origin (tail z))

lookupLocation :: Int -> [(Int, String)] -> String
lookupLocation _ []     = error "No location for that number"
lookupLocation i (x:xs) = if i == fst x then snd x else lookupLocation i xs

locations :: [(Int, String)]
locations = [  (1, "Los Angeles"),
               (2, "San Francisco"),
               (3, "Portland"),
               (4, "Seattle"),
               (5, "Austin"),
               (6, "Trenton"),
               (7, "Nashville"),
               (8, "Lansing"),
               (9, "Denver"),
               (10, "Madison") ]

distanceLookup :: Int -> Int -> [(Int, Int, Int)] -> Int
distanceLookup _ 0 _                     = 0
distanceLookup o d []                    = error ("That origin: "++ show o ++ " /destination "++ show d ++ " pair does not exist")
distanceLookup origin destination (x:xs) = let (f,s,t) = x in if (and [origin == f, destination == s]) then t else distanceLookup origin destination xs

distances :: [(Int, Int, Int)]
distances = [(1, 2, 381),
             (1, 3, 962),
             (1, 4, 1135),
             (1, 5, 1377),
             (1, 6, 2730),
             (1, 7, 2004),
             (1, 8, 2215),
             (1, 9, 1016),
             (1, 10, 1975),
             (3, 1, 962),
	     (3, 2, 635),
             (3, 4, 173),
             (3, 5, 2055),
             (3, 6, 2884),
             (3, 7, 2350),
             (3, 8, 2322),
             (3, 9, 1259),
             (3, 10, 1997),
             (2, 1, 381),
	     (2, 3, 635),
             (2, 4, 807),
             (2, 5, 1757),
             (2, 6, 2895),
             (2, 7, 2303),
             (2, 8, 2332),
             (2, 9, 1269),
             (2, 10, 2092),
             (4, 1, 1135),
	     (4, 3, 173),
             (4, 2, 807),
             (4, 5, 2129),
             (4, 6, 2841),
             (4, 7, 2445),
             (4, 8, 2280),
             (4, 9, 1333),
             (4, 10, 1924),
             (5, 1, 1377),
	     (5, 3, 2055),
             (5, 4, 2129),
             (5, 2, 1757),
             (5, 6, 1700),
             (5, 7, 858),
             (5, 8, 1309),
             (5, 9, 926),
             (5, 10, 1223),
             (6, 1, 2730),
	     (6, 3, 2884),
             (6, 4, 2841),
             (6, 2, 2895),
             (6, 5, 1700),
             (6, 7, 842),
             (6, 8, 668),
             (6, 9, 1749),
             (6, 10, 926),
             (7, 1, 2004),
	     (7, 3, 2350),
             (7, 4, 2445),
             (7, 2, 2303),
             (7, 5, 1700),
             (7, 6, 842),
             (7, 8, 536),
             (7, 9, 1159),
             (7, 10, 623),
             (8, 1, 2215),
	     (8, 3, 2322),
             (8, 4, 2280),
             (8, 2, 2332),
             (8, 5, 1309),
             (8, 6, 668),
             (8, 7, 536),
             (8, 9, 1205),
             (8, 10, 364),
             (9, 1, 1016),
	     (9, 3, 1259),
             (9, 4, 1333),
             (9, 2, 1269),
             (9, 5, 629),
             (9, 6, 1749),
             (9, 7, 1159),
             (9, 8, 1205),
             (9, 10, 964),
             (10, 1, 1975),
	     (10, 3, 1997),
             (10, 4, 1924),
             (10, 2, 2029),
             (10, 5, 1223),
             (10, 6, 926),
             (10, 7, 623),
             (10, 8, 364),
             (10, 9, 964) ]
























