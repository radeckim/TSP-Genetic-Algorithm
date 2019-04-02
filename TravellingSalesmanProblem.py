import scipy.io
import matplotlib.pyplot as draw
import numpy as np
import random
import itertools
from scipy.spatial import distance
from copy import copy
import time
import cProfile
from operator import attrgetter
from scipy.io import loadmat

##Change the working directory
import os
os.chdir("H:\My Desktop")

##Read data 
data = scipy.io.loadmat('xy.mat')

##Class representing a point on the map - a gene
class City(object):

	#create a city by its coordinates and unique ID
    def __init__(self, coord, key):
        self.coord = coord
        self.key = key

##Class representing a chromosome, each gene is acity		
class Route(object):
 
	#Create a chromosome and order it
    def __init__(self, citiesIdxs, citiesList, newMap):
        
        self.citiesIdxs = citiesIdxs
        self.citiesList = citiesList
        self.newMap = newMap
 
        self.orderedList = [self.citiesList[i] for i in self.citiesIdxs]
 
        self.calculateFit()
 
	#Calculate the fitness function - distance
    def calculateFit(self):
 
        self.fitness = self.newMap.calcDist(self.orderedList)
 
        return self.fitness
		
	#randomly swap two indexes
    def swapMutation(self):
 
        pick = random.sample(range(0, 100), 2)
        self.citiesIdxs[pick[0]], self.citiesIdxs[pick[1]] = self.citiesIdxs[pick[1]], self.citiesIdxs[pick[0]]
 
        self.orderedList = [self.citiesList[i] for i in self.citiesIdxs]
	
	#do the inverse mutation
    def inversionMutation(self):
    
        pick = random.sample(range(0, 100), 2)
    
        if(pick[0] > pick[1]):
    
            temp = pick[1]
    
            pick[1] = pick[0]
            pick[0] = temp
    
        self.orderedList[pick[0]:pick[1]] = self.orderedList[pick[0]:pick[1]][::-1]
	
	#do the flip mutation
    def flipMutation(self):

        pick = random.sample(range(0, 100), 2) 

        if(pick[0] > pick[1]):

            temp = pick[1]
   
            pick[1] = pick[0]
            pick[0] = temp
        
        self.orderedList.insert(pick[1], self.orderedList[pick[0]])
        self.orderedList.pop(pick[0])      
	
	#diplay an order of a particular chromosome
    def dispOrder(self):
 
        return [city.key for city in self.orderedList]
	
	#get fitness of a particular chromosome
    def getFitness(self):
 
        return self.fitness
    
	#display a map with linked cities in the order
    def displayMap(self):
 
        draw.plot([city.coord[0] for city in self.orderedList],[city.coord[1] for city in self.orderedList], 'xb-')
        draw.show()

##a class representing a whole area on which the cities are		
class Map(object):

	#initilialize the map by calculating a distances between each two cities 
    def __init__(self, citiesList):
 
        self.citiesList = citiesList
        self.cityDistances = {}
 
        for x in self.citiesList:
            for y in self.citiesList:
                if x == y:
                    continue
                self.cityDistances[(x, y)]= ((x.coord[0] - y.coord[0]) ** 2 + (x.coord[1] - y.coord[1]) ** 2) ** 0.5 
	
	#all of the distances can be reached from that function - there's no need to recalculate them
    def calcDist(self, orderedCities):
 
        self.distance = 0
 
        for i, n in enumerate(orderedCities):
 
            if(i == len(orderedCities) - 1):
 
                self.distance = self.distance + self.cityDistances[(n, orderedCities[0])]
 
            else:
 
                self.distance = self.distance + self.cityDistances[(n, orderedCities[i + 1])]
 
        return self.distance

##a class representing the whole flow of the GA 		
class Salesman(object):
 
	##define the population size
    def __init__(self, populationSize):
 
        self.populationSize = populationSize        
	
	##load data from given file
    def loadData(self, data):
 
        rawList = data.get('xy').tolist();
        self.citiesList = [City(row, rawList.index(row)) for row in rawList]
 
        self.newMap = Map(self.citiesList)
	
	##initial random population 
    def populate(self):
 
        self.population = [Route(random.sample(range(0, 100), 100), self.citiesList, self.newMap) for x in
                           range(self.populationSize)]
	
	##get fitness of the whole population
    def getPopulation(self):
 
        return [c.getFitness() for c in self.population]
	
	##get the fittest individual of the current population
    def getFittest(self):

        return min(self.population, key=attrgetter('fitness'))
	
	##get list of ordered cities
    def getCitiesList(self):
 
        return self.citiesList
	
	##roulette selection implementation
    def selectRoulette(self):
 
        max = sum(self.getPopulation())
        pick = random.uniform(0, max)
 
        current = 0
 
        for chromosome in self.population:
 
            current += chromosome.getFitness()
            if current > pick:
                return chromosome
				
	##tournament selection implementation
    def selectTournament(self, iterations, size):
 
        aspirants = [random.choice(self.population) for i in range(size)]
        chosen = (min(aspirants, key=attrgetter('fitness')))
        return chosen
	
	##select two parents
    def selectTwo(self):
 
        while(True):
 
            self.road1 = self.selectTournament(1, 4)
            self.road2 = self.selectTournament(1, 4)
			
			#make sure that the parents are not the same
            if(self.road1 == self.road2):
 
                continue
 
            else:
 
                break 
	##order one permuation 
    def orderOnePermut(self, winner, looser):
 
        self.winnerOrder = winner
        self.looserOrder = looser
 
        a, b = random.sample(range(len(looser)), 2)
		
		#make sure the first number is smaller then the first one
        if a > b:
 
            temp = b
            b = a
            a = temp
 
        self.child1 = [None] * 100
		
		#choose a substring
        self.child1[a:b] = self.winnerOrder[a:b]
 
		#find all of the elements which are in the parent2 but not in the substring
        uniques1 = [item for item in self.looserOrder if item not in self.child1] 
        position = 0
		
		#place all of the missing values	
        for i, n in enumerate(self.child1):
 
            if n is None:
 
                self.child1[i] = uniques1[position]
                position = position + 1
 
        return self.child1
	
	#the actual flow of the algorithm
    def evaluatePopulation(self, eliteSize):

        #define constance values and declare new population array    
        crossoverConst = 0.5
        mutateConst = 0.05
        newPopulation = []

        #evaluate n times (n = populationSize) 
        for i in range(self.populationSize):

            #select parents    
            self.selectTwo()
 
            first = self.road1.dispOrder()
            second = self.road2.dispOrder()

            #decide if crossover occurs
            if (random.random() < crossoverConst):
 
                newChromosome1 = self.orderOnePermut(first, second)
 
            else:

                #if not don't change parents 
                newChromosome1 = first

            #append a new member
            newPopulation.append(Route(newChromosome1, self.citiesList, self.newMap))

        #random mutations - 3 types in a row
        for i in range(self.populationSize):
 
            if(random.random() < mutateConst):
 
                newPopulation[i].swapMutation()

            if random.random() < mutateConst:                
   
               newPopulation[i].inversionMutation()      

            if random.random() < mutateConst:              
    
               newPopulation[i].flipMutation()     

        #override the old population
        self.population = newPopulation 

def main():

    #Test of the algorithm 
    x = Salesman(200) 
    x.loadData(data) 
    x.populate()
 
    i = 0
    b = []

    #Run the algorithm for 1k iterations 
    while i < 10000:
 
        x.evaluatePopulation(0)
        b.append(min(x.getPopulation()))
        i = i + 1
        print(i, min(x.getPopulation()))     
  
    #plot a graph with each generation and its minimal fitness score    
    draw.plot(range(i), b)
    draw.show()

    #show the chromosome with the minimal fitness score in the last generation
    print(min(b))   
 
main()
