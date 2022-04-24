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


SOLVERS: Dict[str, Callable[[Instance], Solution]] = {
    "naive": solve_naive, "greedySetCover": solve_greedySetCover
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
            #print("# Penalty: ", solution.penalty(), file=g)
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
