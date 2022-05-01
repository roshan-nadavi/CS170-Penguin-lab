"""Solves an instance.

Modify this file to implement your own solvers.

For usage, run `python3 solve.py --help`.
"""

import argparse
from pathlib import Path
from typing import Callable, Dict

from instance import Instance
from solution import Solution
from file_wrappers import StdinFileWrapper, StdoutFileWrapper
from point import Point
#from point import distance_sq


def solve_naive(instance: Instance) -> Solution:
    return Solution(
        instance=instance,
        towers=instance.cities,
    )

def solve_greedySetCoverWCaveat(instance: Instance) -> Solution:
    #initialize the sets
    pointCoverage = []
    if instance.grid_side_length == 30:
        pointCoverage = [[] for x in range(900)]
    if instance.grid_side_length == 50:
        pointCoverage = [[] for x in range(50*50)]
    if instance.grid_side_length == 100:
        pointCoverage = [[] for x in range(100*100)]
    for x in range(instance.grid_side_length):
        for y in range(instance.grid_side_length):
            for i in range(x-instance.R_s, x+instance.R_s):
                for j in range(y-instance.R_s, y+instance.R_s):
                    if ((Point(x, y).distance_sq(Point(i,j))**0.5 )) <= instance.R_s:
                        if Point(i, j) in instance.cities:
                            pointCoverage[x+y*instance.grid_side_length].append(Point(i,j))
    citiesLeft = instance.cities.copy()
    pointsUsed = []
        
    while len(citiesLeft)>0:
        #print(len(citiesLeft))
        optPoint = 0
        coverage = len(pointCoverage[0])
        for i in range(len(pointCoverage)):
            if len(pointCoverage[i]) >= coverage:
                if len(pointCoverage[i]) >= coverage or i < instance.grid_side_length **2:
                    optPoint = i
                    coverage = len(pointCoverage[i])
        pointsUsed.append(Point(optPoint%instance.grid_side_length, 
                                optPoint//instance.grid_side_length))
        #now remove all the points that are covered
        temp = pointCoverage[optPoint].copy()
        #print(temp, optPoint)
        for singPoint in temp:
            for singleList in pointCoverage:
                if singPoint in singleList:
                    singleList.remove(singPoint)
                if singPoint in citiesLeft:
                    citiesLeft.remove(singPoint)   
    #return solution
    return Solution(
        instance=instance,
        towers=pointsUsed,
    )

def solve_greedySetCoverWOPenalty(instance: Instance) -> Solution:
    #initialize the sets
    pointCoverage = []
    if instance.grid_side_length == 30:
        pointCoverage = [[] for x in range(900)]
    if instance.grid_side_length == 50:
        pointCoverage = [[] for x in range(50*50)]
    if instance.grid_side_length == 100:
        pointCoverage = [[] for x in range(100*100)]
    for x in range(instance.grid_side_length):
        for y in range(instance.grid_side_length):
            for i in range(x-instance.R_s, x+instance.R_s):
                for j in range(y-instance.R_s, y+instance.R_s):
                    if ((Point(x, y).distance_sq(Point(i,j))**0.5 )) <= instance.R_s:
                        if Point(i, j) in instance.cities:
                            pointCoverage[x+y*instance.grid_side_length].append(Point(i,j))
    citiesLeft = instance.cities.copy()
    pointsUsed = []
        
    while len(citiesLeft)>0:
        #print(len(citiesLeft))
        optPoint = 0
        coverage = len(pointCoverage[0])
        for i in range(len(pointCoverage)):
            if len(pointCoverage[i]) > coverage:
                optPoint = i
                coverage = len(pointCoverage[i])
        pointsUsed.append(Point(optPoint%instance.grid_side_length, 
                                optPoint//instance.grid_side_length))
        #now remove all the points that are covered
        temp = pointCoverage[optPoint].copy()
        #print(temp, optPoint)
        for singPoint in temp:
            for singleList in pointCoverage:
                if singPoint in singleList:
                    singleList.remove(singPoint)
                if singPoint in citiesLeft:
                    citiesLeft.remove(singPoint)   
    #return solution
    return Solution(
        instance=instance,
        towers=pointsUsed,
    )
def solve_greedySetCover(instance: Instance) -> Solution:
    #initialize the sets
    pointCoverage = []
    if instance.grid_side_length == 30:
        pointCoverage = [[] for x in range(900)]
    if instance.grid_side_length == 50:
        pointCoverage = [[] for x in range(50*50)]
    if instance.grid_side_length == 100:
        pointCoverage = [[] for x in range(100*100)]
    for x in range(instance.grid_side_length):
        for y in range(instance.grid_side_length):
            for i in range(x-instance.R_s, x+instance.R_s):
                for j in range(y-instance.R_s, y+instance.R_s):
                    if ((Point(x, y).distance_sq(Point(i,j))**0.5 )) <= instance.R_s:
                        if Point(i, j) in instance.cities:
                            pointCoverage[x+y*instance.grid_side_length].append(Point(i,j))
    citiesLeft = instance.cities.copy()
    pointsUsed = []
    #midpoint = (instance.grid_side_length*instance.grid_side_length)/2
    def overlaps(currPoint):
        overlap = 0
        for aPoint in pointsUsed:
            if (Point(currPoint%instance.grid_side_length, 
                                currPoint//instance.grid_side_length).distance_sq(aPoint)) ** 0.5 <= instance.R_p:
                overlap = overlap + 1;
        return overlap
            
            
    while len(citiesLeft)>0:
        #print("iteration")
        #print(len(citiesLeft))
        optPoint = 0
        coverage = len(pointCoverage[0])
        penalty = overlaps(0)
        for i in range(len(pointCoverage)):
            pointCovLength = len(pointCoverage[i])
            if pointCovLength >= coverage:
                tempPenalty = overlaps(i)
                if tempPenalty < penalty or pointCovLength > coverage:
                    optPoint = i
                    coverage = len(pointCoverage[i])
                    penalty = tempPenalty
        pointsUsed.append(Point(optPoint%instance.grid_side_length, 
                                optPoint//instance.grid_side_length))
        #now remove all the points that are covered
        temp = pointCoverage[optPoint].copy()
        #print(temp, optPoint)
        for singPoint in temp:
            for singleList in pointCoverage:
                if singPoint in singleList:
                    singleList.remove(singPoint)
            if singPoint in citiesLeft:
                citiesLeft.remove(singPoint)   
    #return solution
    return Solution(
        instance=instance,
        towers=pointsUsed,
    )

def solve_greedyTwoSetCover(instance: Instance) -> Solution:
    #initialize the sets
    dim = instance.grid_side_length
    pointCoverage = []
    if instance.grid_side_length == 30:
        pointCoverage = [[] for x in range(900)]
    if instance.grid_side_length == 50:
        pointCoverage = [[] for x in range(50*50)]
    if instance.grid_side_length == 100:
        pointCoverage = [[] for x in range(100*100)]
    for x in range(instance.grid_side_length):
        for y in range(instance.grid_side_length):
            for i in range(x-instance.R_s, x+instance.R_s):
                for j in range(y-instance.R_s, y+instance.R_s):
                    if ((Point(x, y).distance_sq(Point(i,j))**0.5 )) <= instance.R_s:
                        if Point(i, j) in instance.cities:
                            pointCoverage[x+y*instance.grid_side_length].append(Point(i,j))
    citiesLeft = instance.cities.copy()
    pointsUsed = []
    #midpoint = (instance.grid_side_length*instance.grid_side_length)/2
    def overlaps(currPoint, otherPoint):
        overlap = 0
        if currPoint == otherPoint:
            overlap = overlap - 0.5
            for aPoint in pointsUsed:
                if (Point(currPoint%instance.grid_side_length, 
                                    currPoint//instance.grid_side_length).distance_sq(aPoint)) ** 0.5 <= instance.R_p:
                    overlap = overlap + 1;
        else:
             for aPoint in pointsUsed:
                if (Point(currPoint%instance.grid_side_length, 
                                    currPoint//instance.grid_side_length).distance_sq(aPoint)) ** 0.5 <= instance.R_p:
                    overlap = overlap + 1;
                if (Point(otherPoint%instance.grid_side_length, 
                                    otherPoint//instance.grid_side_length).distance_sq(aPoint)) ** 0.5 <= instance.R_p:
                    overlap = overlap + 1;
             aPoint = Point(otherPoint%dim, otherPoint//dim)
             if (Point(currPoint%dim, currPoint//dim).distance_sq(aPoint)) ** 0.5 <= instance.R_p:
                    overlap = overlap + 1
        return overlap
            
            
    while len(citiesLeft)>0:
        #print("iteration")
        #print(len(citiesLeft))
        optPoint = 0
        optPoint2 = 0
        coverage = len(pointCoverage[0])
        penalty = overlaps(0,0)
        for i in range(len(pointCoverage)):
            for j in range(i, len(pointCoverage)):
                #print(i, j)
                multilist = pointCoverage[i] + pointCoverage[j]
                multiset = set(multilist)
                multilist = list(multiset)
                pointCovLength = len(multilist)
                if pointCovLength >= coverage:
                    tempPenalty = overlaps(i, j)
                    if tempPenalty < penalty or pointCovLength > coverage:
                        optPoint = i
                        optPoint2 = j
                        coverage = len(multilist)
                        penalty = tempPenalty
        pointsUsed.append(Point(optPoint%instance.grid_side_length, 
                                optPoint//instance.grid_side_length))
        pointsUsed.append(Point(optPoint2%instance.grid_side_length, 
                                optPoint2//instance.grid_side_length))
        #now remove all the points that are covered
        temp = pointCoverage[optPoint].copy()
        temp2 = pointCoverage[optPoint2].copy()
        #print(temp, optPoint)
        for singPoint in temp:
            for singleList in pointCoverage:
                if singPoint in singleList:
                    singleList.remove(singPoint)
            if singPoint in citiesLeft:
                citiesLeft.remove(singPoint) 
        for singPoint in temp2:
            for singleList in pointCoverage:
                if singPoint in singleList:
                    singleList.remove(singPoint)
            if singPoint in citiesLeft:
                citiesLeft.remove(singPoint) 
    #return solution
    pointsUsedSet = set(pointsUsed)
    pointsUsed = list(pointsUsedSet)
    return Solution(
        instance=instance,
        towers=pointsUsed,
    )
def solve_naive_better(instance: Instance) -> Solution:
    tows = instance.cities
    for i in range(len(tows)):
        for j in range(len(tows)):
            if tows[i].distance_sq(tows[j])**0.5 <= instance.R_s and i != j:
                tows[j] = Point(-100, -100)
    while Point(-100, -100) in tows:
        tows.remove(Point(-100, -100))
    return Solution(
        instance=instance,
        towers=tows,
    )
SOLVERS: Dict[str, Callable[[Instance], Solution]] = {
    "naive": solve_naive, "greedySetCover": solve_greedySetCover,
    "greedyTwoSetCover": solve_greedyTwoSetCover,
    "greedySetCoverWOPenalty":solve_greedySetCoverWOPenalty,
    "greedySetCoverWCaveat":solve_greedySetCoverWCaveat,
    "betterNaive":solve_naive_better
}


# You shouldn't need to modify anything below this line.
def infile(args):
    if args.input == "-":
        return StdinFileWrapper()

    return Path(args.input).open("r")


def outfile(args):
    if args.output == "-":
        return StdoutFileWrapper()

    return Path(args.output).open("w")


def main(args):
    with infile(args) as f:
        instance = Instance.parse(f.readlines())
        solver = SOLVERS[args.solver]
        solution = solver(instance)
        assert solution.valid()
        with outfile(args) as g:
            print("# Penalty: ", solution.penalty(), file=g)
            solution.serialize(g)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Solve a problem instance.")
    parser.add_argument("input", type=str, help="The input instance file to "
                        "read an instance from. Use - for stdin.")
    parser.add_argument("--solver", required=True, type=str,
                        help="The solver type.", choices=SOLVERS.keys())
    parser.add_argument("output", type=str,
                        help="The output file. Use - for stdout.",
                        default="-")
    main(parser.parse_args())
