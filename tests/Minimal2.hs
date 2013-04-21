import Genetics2
import qualified Data.ByteString as B

config = defaultConfig { geneSize = 4,
                         chromosomeSize = 5,
                         populationSize = 100,
                         targetFitness  = 100,
                         crossoverRate  = 0.9,
                         mutationRate   = 0.001,
                         fitnessFunction = decode,
                         maxGen           = 100
                        }

decode :: BitString -> Float
decode (BitString bs) = B.foldl (\acc -> \v -> if (v == 1) then acc+5 else acc) 0 bs

