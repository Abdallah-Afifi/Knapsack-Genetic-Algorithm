#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <iostream>
#include <algorithm>

using namespace std;

#define RANDOM_DECIMAL (rand() / (double) RAND_MAX)

#define POPULATION_SIZE 200
#define NUM_GENERATIONS 400
#define ELITE_COUNT 1
const double CROSSOVER_PROBABILITY = 0.85;
const double MUTATION_PROBABILITY = 0.1;

struct Chromosome
{
    unsigned int genes;
    int fitness;
    int weight;
    Chromosome() {}
};

int numGenes, knapsackSize;
int weights[55];
int values[55];
double cumulativeDistributionFunction[POPULATION_SIZE];

Chromosome currentPopulation[POPULATION_SIZE];
Chromosome nextPopulation[POPULATION_SIZE];

void initialize()
{
    memset(currentPopulation, 0, sizeof(currentPopulation));
    srand(time(0));
    for (int i = 0; i < POPULATION_SIZE; i++)
    {
        int currentWeight = 0;
        for (int j = 0; j < numGenes && currentWeight < knapsackSize; j++)
        {
            if (rand() % 2 && currentWeight + weights[j] <= knapsackSize)
            {
                currentPopulation[i].genes |= (1 << j);
                currentWeight += weights[j];
            }
        }
    }
}

void calculateFitness()
{
    int totalFitness = 0;
    memset(cumulativeDistributionFunction, 0, sizeof(cumulativeDistributionFunction));
    for (int i = 0; i < POPULATION_SIZE; i++)
    {
        currentPopulation[i].fitness = 0;
        currentPopulation[i].weight = 0;
        for (int j = 0; j < numGenes; j++)
        {
            currentPopulation[i].fitness += ((currentPopulation[i].genes >> j) & 1) * values[j];
            currentPopulation[i].weight += ((currentPopulation[i].genes >> j) & 1) * weights[j];
        }
        if (currentPopulation[i].weight > knapsackSize)
            currentPopulation[i].fitness = 0;
        totalFitness += currentPopulation[i].fitness;
    }
    for (int i = 0; i < POPULATION_SIZE; i++)
    {
        double probabilityOfI = currentPopulation[i].fitness / (double)totalFitness;
        cumulativeDistributionFunction[i] = (i > 0 ? cumulativeDistributionFunction[i - 1] + probabilityOfI : probabilityOfI);
    }
}

bool isValid(Chromosome *pChromosome)
{
    int totalWeight = 0;
    for (int i = 0; i < numGenes; i++)
        totalWeight += ((pChromosome->genes >> i) & 1) * weights[i];
    return totalWeight <= knapsackSize;
}

void evolve()
{
    int currentIndex = 0;
    if (ELITE_COUNT)
    {
        pair<int, int> best(0, 0), secondBest(0, 0);
        for (int i = 0; i < POPULATION_SIZE; i++)
        {
            if (currentPopulation[i].fitness >= best.second)
            {
                secondBest = best;
                best.first = i;
                best.second = currentPopulation[i].fitness;
            }
            else if (currentPopulation[i].fitness > secondBest.second)
            {
                secondBest.first = i;
                secondBest.second = currentPopulation[i].fitness;
            }
        }
        nextPopulation[currentIndex++] = currentPopulation[best.first];
        nextPopulation[currentIndex++] = currentPopulation[secondBest.first];
    }

    while (currentIndex < POPULATION_SIZE)
    {
        Chromosome parent1 = currentPopulation[upper_bound(cumulativeDistributionFunction, cumulativeDistributionFunction + POPULATION_SIZE, RANDOM_DECIMAL) - cumulativeDistributionFunction];
        Chromosome parent2 = currentPopulation[upper_bound(cumulativeDistributionFunction, cumulativeDistributionFunction + POPULATION_SIZE, RANDOM_DECIMAL) - cumulativeDistributionFunction];
        Chromosome child1 = parent1;
        Chromosome child2 = parent2;
        nextPopulation[currentIndex++] = parent1;
        nextPopulation[currentIndex++] = parent2;
        bool willCrossover = (RANDOM_DECIMAL <= CROSSOVER_PROBABILITY);
        if (willCrossover)
        {
            int crossoverPoint = rand() % numGenes;
            if (crossoverPoint == 0)
                crossoverPoint++;
            else if (crossoverPoint == numGenes)
                crossoverPoint--;
            Chromosome child1 = parent1;
            Chromosome child2 = parent2;
            for (int i = crossoverPoint + 1; i < numGenes; i++)
            {
                if (parent1.genes & (1 << (i - 1)))
                    child2.genes |= 1 << (i - 1);
                else
                    child2.genes &= ~(1 << (i - 1));

                if (parent2.genes & (1 << (i - 1)))
                    child1.genes |= 1 << (i - 1);
                else
                    child1.genes &= ~(1 << (i - 1));
            }
        }
        if (isValid(&child1) && currentIndex < POPULATION_SIZE)
            nextPopulation[currentIndex++] = child1;
        if (isValid(&child2) && currentIndex < POPULATION_SIZE)
            nextPopulation[currentIndex++] = child2;
    }

    for (int i = 0; i < POPULATION_SIZE; i++)
    {
        if (RANDOM_DECIMAL <= MUTATION_PROBABILITY)
        {
            int flip = rand() % numGenes;
            Chromosome mutated = nextPopulation[i];
            mutated.genes ^= (1 << flip);
            if (isValid(&mutated))
                nextPopulation[i] = mutated;
        }
    }
}

int main()
{
    int cases;
    cin >> cases;
    while (cases--)
    {
        cin >> numGenes >> knapsackSize;
        for (int i = 0; i < numGenes; i++)
            cin >> weights[i] >> values[i];
        initialize();
        for (int i = 0; i < NUM_GENERATIONS; i++)
        {
            calculateFitness();
            evolve();
            memcpy(currentPopulation, nextPopulation, sizeof(nextPopulation));
        }
        int globalOptimum = 0, globalOptimumIndex = 0;
        calculateFitness();
        for (int i = 0; i < POPULATION_SIZE; i++)
        {
            if (currentPopulation[i].fitness > globalOptimum)
            {
                globalOptimum = currentPopulation[i].fitness;
                globalOptimumIndex = i;
            }
        }
        for (int i = 0; i < numGenes; i++)
            cout << ((currentPopulation[globalOptimumIndex].genes >> i) & 1);
        cout << " = " << globalOptimum << endl;
    }
    return 0;
}
