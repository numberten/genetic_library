Genetic Library
===============


The idea behind this project is to create a library that hides the algorithm aspect of genetic algorithms. Using a genetic algorithm to solve a problem should be as easy as putting the problem in terms of a genetic algorithm (which in itself isn't always an easy task). With this module genetic algorithms should be as simple as declaring a configuration struct with appropriately assigned values and passing it to the runGA function. 

Currently Configurable Values
-----------------------------

   - geneSize        (the size, in bits, of a gene)
   - chromosomeSize  (how many genes each chromosome will have)
   - populationSize  (the number of chromosomes your population will consist of)
   - maxGen          (a limit as to how many times generations the algorithm will breed)
   - targetFitness   (the fitness goal upon which the algorithm ends)
   - crossoverRate   (the likelihood that genes selected to breed will crossover)
   - fitnessFunction (a function that takes a chromosome and returns a fitness)
   - ordered         (boolean value, toggle for ordered chromosomes)
   - mutationRate    (the probability that a gene will be selected for mutation)

Configurable Values To Come
---------------------------

   - ~~mutationRate~~
   - flag for two-point crossover
   - ~~flag for ordered chromosomes~~
   - tournament selection option
   - option for user defined genes
