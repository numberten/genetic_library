import Genetics

config = defaultConfig { geneSize = 4,
                         chromosomeSize = 5,
                         populationSize = 100,
                         targetFitness  = 100,
                         crossoverRate  = 0.9,
                         mutationRate   = 0.001,
                         fitnessFunction = decode,
                         maxGen           = 100
                        }

decode :: Chromosome -> Float
decode bs = foldl (\acc -> \v -> acc + (accum 0 v)) 0 bs

accum :: Float -> Gene -> Float
accum f (BitString xs) = foldl (\acc -> \v -> if (v == '1') then acc+5 else acc) f xs

