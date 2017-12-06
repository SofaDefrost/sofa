import geneticalgo
import random
import copy
import os
import json
import sys
import math
import crochetutils
from geneticalgo import algorithm
from geneticalgo import primitives
from geneticalgo import primitives2D
from geneticalgo import piecewisePolynom2D
from geneticalgo import crochetInFile
from sofalauncher import launcher
import crochet
import time
import errno


A = [0.0, 0.010, ["smooth", 0.0007, 0.0], ["smooth", 0.0007, 0.0]]
B = [0.0, 0.005, ["smooth", 0.0007, 1.0], ["smooth", 0.0007, 1.0]]
C = [-0.0025, 0.00125, ["smooth", 0.0007, 1.0], ["smooth", 0.0007, 1.0]]
D = [-0.0025, -0.00125, ["smooth", 0.0007, 1.0], ["corner", 0.0015, 1.0, 0.0007, 0.0007]]
E = [0.0, -0.00250, ["smooth", 0.0007, 1.0], ["smooth", 0.0007, 1.0]]
F = [0.0025, -0.0025, ["smooth", 0.0007, 1.0], ["smooth", 0.0007, 1.0]]
G = [0.0035, -0.00125, ["smooth", 0.0007, 1.0], ["smooth", 0.0007, 1.0]]

reference = [A, B, C, D, E, F, G]

individualId = 0
mutationType="ON"
mutationThickness="ON"
mutationCoef="OFF"
mutationWidth="OFF"
mutationDepht="OFF"
mutationPosition="ON"

def getNextId():
    global individualId
    individualId+=1
    return individualId

def getJSONFragmentFrom(file):

    with open(file, "r") as source:
        for line in source:
            if "JSON : " in line:
                data = json.loads(line[7:])
                return data

    return None




class CrochetIndividual(algorithm.Individual):

    def __init__(self):
        self.level=None
        self.id = getNextId()
        crochetutils.create(self)



def newIndividualFrom(ind):
    newInd=CrochetIndividual()
    newInd.listOfDrawnPoints=copy.deepcopy(ind.listOfDrawnPoints)

    return newInd

def getShapeFromInd(ind):
    return crochet.generateCrochetManually(ind.listOfDrawnPoints)



####################################################################################################
### MUTATION
###



def mutation_Position(ind, side):

    length=len(ind.listOfDrawnPoints)
    if length < 2:
        raise ValueError, "don't touch the bottom and the top!"

    if not side in ["left", "right"]:
        raise ValueError, "choose a side"
    index=random.randint(1,length-2)
    ind.listOfDrawnPoints[index][0]+=random.uniform(-0.0001,0.0001)
    ind.listOfDrawnPoints[index][1]+=random.uniform(-0.0001,0.0001)





def mutation_Depht(ind, side):

    length=len(ind.listOfDrawnPoints)
    if length < 2:
        raise ValueError, "don't touch the bottom and the top!"

    if not side in ["left", "right"]:
        raise ValueError, "choose a side"


    index=random.randint(1,length-2)

    if side == "left":

        param = ind.listOfDrawnPoints[index][2]

        if param[0]=="corner":
            Depht=param[4]
            epsilon = random.uniform(0.8,1.2)
            ind.listOfDrawnPoints[index][2][4]=min(ind.listOfDrawnPoints[index][2][4]*epsilon, param[1]/3.0)
#            ind.listOfDrawnPoints[index][2][2]*=1.0/epsilon
    else:

        param = ind.listOfDrawnPoints[index][3]

        if param[0]=="corner":
            Depht=param[4]
            epsilon = random.uniform(0.8,1.2)
            ind.listOfDrawnPoints[index][3][4]=min(ind.listOfDrawnPoints[index][2][4]*epsilon, param[1]/3.0)
#            ind.listOfDrawnPoints[index][3][2]*=1.0/epsilon

def mutation_Width(ind, side):

    length=len(ind.listOfDrawnPoints)
    if length < 2:
        raise ValueError, "don't touch the bottom and the top!"

    if not side in ["left", "right"]:
        raise ValueError, "choose a side"


    index=random.randint(1,length-2)

    if side == "left":

        param = ind.listOfDrawnPoints[index][2]

        if param[0]=="corner":
            width=param[3]
            epsilon = random.uniform(0.8,1.2)
            ind.listOfDrawnPoints[index][2][3]*=epsilon

    else:

        param = ind.listOfDrawnPoints[index][3]

        if param[0]=="corner":
            width=param[3]
            epsilon = random.uniform(0.8,1.2)
            ind.listOfDrawnPoints[index][3][3]*=epsilon

def mutation_Type(ind, side):

    length=len(ind.listOfDrawnPoints)
    if length < 2:
        raise ValueError, "don't touch the bottom and the top!"

    if not side in ["left", "right"]:
        raise ValueError, "choose a side"

    index=random.randint(1,length-2)

    if side == "left":

        param = ind.listOfDrawnPoints[index][2]
        type=param[0]

        if type == "smooth":

            type = "corner"
            ind.listOfDrawnPoints[index][2][0] = type
            ind.listOfDrawnPoints[index][2].append(min(param[1]/3.0,0.0003))
            ind.listOfDrawnPoints[index][2].append(min(param[1]/3.0,0.0003))
            ind.listOfDrawnPoints[index][2][2]*=0.5

            print str(ind.listOfDrawnPoints[index][2])
        else:
            type = "smooth"
            ind.listOfDrawnPoints[index][2][0] = type
            del ind.listOfDrawnPoints[index][2][-2:]
            print str(ind.listOfDrawnPoints[index][2])

    else:

        param = ind.listOfDrawnPoints[index][3]
        type=param[0]

        if type == "smooth":

            type = "corner"
            ind.listOfDrawnPoints[index][3][0] = type
            ind.listOfDrawnPoints[index][3].append(min(param[1]/3.0,0.0003))
            ind.listOfDrawnPoints[index][3].append(min(param[1]/3.0,0.0003))
            ind.listOfDrawnPoints[index][3][2]*=0.5
            print str(ind.listOfDrawnPoints[index][3])
        else:
            type = "smooth"
            ind.listOfDrawnPoints[index][3][0] = type
            del ind.listOfDrawnPoints[index][3][-2:]
            print str(ind.listOfDrawnPoints[index][3])


def mutation_Thickness(ind, side):

    length=len(ind.listOfDrawnPoints)

    if not side in ["left", "right"]:
        raise ValueError, "choose a side"

    index=random.randint(0,length-1)

    if side == "left":

        param = ind.listOfDrawnPoints[index][2]
        thickness=param[1]
        epsilon=random.uniform(0.9, 1.1)
        thickness*=epsilon
        ind.listOfDrawnPoints[index][2][1] = thickness
    else:
        param = ind.listOfDrawnPoints[index][3]
        thickness=param[1]
        epsilon=random.uniform(0.9, 1.1)
        thickness*=epsilon
        ind.listOfDrawnPoints[index][3][1] = thickness



def mutation_Coef(ind, side):

    length=len(ind.listOfDrawnPoints)
    if length < 2:
        raise ValueError, "don't touch the bottom and the top!"

    if not side in ["left", "right"]:
        raise ValueError, "choose a side"

    index=random.randint(1,length-2)
    if side == "left":

        param = ind.listOfDrawnPoints[index][2]
        coef = param[2]
        epsilon=random.uniform(0.8,1.2)
        coef=max(0.0, min(1.1, coef*epsilon))
        ind.listOfDrawnPoints[index][2][2] = coef

    else:

        param = ind.listOfDrawnPoints[index][3]
        coef = param[2]
        epsilon=random.uniform(0.8,1.2)
        coef=max(0.0, min(1.1, coef*epsilon))
        ind.listOfDrawnPoints[index][3][2] = coef

def mutation(ind):


    if mutationType=="ON":
        if random.choice([True, False]):
            side = random.choice(["left", "right"])
            mutation_Type(ind, side)

    if mutationThickness=="ON":
        if random.choice([True, False]):
            side = random.choice(["left", "right"])
            mutation_Thickness(ind, side)


    if mutationCoef=="ON":
        if random.choice([True, False]):
            side = random.choice(["left", "right"])
            mutation_Coef(ind, side)

    if mutationWidth=="ON":
        if random.choice([True, False]):
            side = random.choice(["left", "right"])
            mutation_Width(ind, side)

    if mutationDepht=="ON":
        if random.choice([True, False]):
            side = random.choice(["left", "right"])
            mutation_Depht(ind, side)


    if mutationPosition=="ON":
        if random.choice([True, False]):
            side = random.choice(["left", "right"])
            mutation_Position=(ind, side)


def mutationFunc(pop, params):

    number_of_mutated_ind=params["nbIndMutated"]
    number_of_mutation_per_ind=params["nbMutationsPerInd"]

    if not (isinstance(number_of_mutated_ind, int) and isinstance(number_of_mutation_per_ind, int)):
        raise ValueError, "nbIndMutated and nbMutationsPerInd are not given or are not int"

    length_temp=len(pop)

    for i in range(number_of_mutated_ind):
        j=random.randint(0,length_temp-1)
        ind=newIndividualFrom(pop[j])

        for k in range(number_of_mutation_per_ind):
            mutation(ind)

        pop.append(ind)

    return pop

####################################################################################################
###CROSSING
###
def crossing_ind(individual1, individual2):


    crochet1 = crochet.generateCrochetManually(individual1.listOfDrawnPoints)
    crochet2 = crochet.generateCrochetManually(individual2.listOfDrawnPoints)

    s = random.uniform(0.0, 1.0)

    index1 = 0
    while crochet1.listOfControlPoints[index1].absisse < s:
        index1+=1

    index2 = 0
    while crochet2.listOfControlPoints[index2].absisse < s:
        index2+=1


    ind1=newIndividualFrom(individual1)
    ind2=newIndividualFrom(individual2)

    del ind1.listOfDrawnPoints[index1:]
    del ind2.listOfDrawnPoints[index2:]

    ind1.listOfDrawnPoints+=individual2.listOfDrawnPoints[index2:]
    ind2.listOfDrawnPoints+=individual1.listOfDrawnPoints[index1:]

    return (ind1, ind2)


def crossFunc(pop, params):

    number_of_crossing=params["crossTx"]

    if not isinstance(number_of_crossing, int):
        raise ValueError, "crossTx is not given or is not an int"
    print("crossFunc "+str(pop))

    length_temp=len(pop)

    newpop = algorithm.Population()
    for i in range(number_of_crossing):
        j=random.randint(0, length_temp-1)
        k=random.randint(0, length_temp-1)
        (a,b)=crossing_ind(pop[j], pop[k])
        newpop.append(a)
        newpop.append(b)
    return newpop

####################################################################################################
### Generate
###
def generateIndividual(reference):
        individual=CrochetIndividual()
        individual.listOfDrawnPoints = reference

#        side = random.choice(["left", "right"])
#        mutation_Type(individual, side)

        for i in range(5):
            side = random.choice(["left", "right"])
            mutation_Thickness(individual, side)


        for i in range(5):
            side = random.choice(["left", "right"])
            mutation_Coef(individual, side)

#            if mutationWidth=="ON":
#                if random.choice([True, False]):
#                    side = random.choice(["left", "right"])
#                    mutation_Width(individual, side)

#            if mutationDepht=="ON":
#                if random.choice([True, False]):
#                    side = random.choice(["left", "right"])
#                    mutation_Depht(individual, side)


        for i in range(5):
            side = random.choice(["left", "right"])
            mutation_Position=(individual, side)

        return individual
#        length = len(reference)

#        individual.listOfDrawnPoints.append(reference[0])

#        if length > 2:

#            for i in range(1, length-1):

#                param = reference[i]

#                newPoint = [param[0]+random.uniform(-0.0001, 0.0001), param[1]+random.uniform(-0.0001, 0.0001)]

#                if param[2][0] == "smooth":
#                    newPoint.append(["smooth", max(0.0001, param[2][1]+random.uniform(-0.0001, 0.0001), max(0.0, min(1.1,param[2][2]+random.uniform(-0.0001, 0.00001)))])

#                elif param[2][0] == "corner":
#                    newPoint.append(["corner", max(0.0001, param[2][1]+random.uniform(-0.0001, 0.0001)), max(0.0, min(1.1,param[2][2]+random.uniform(-0.0001, 0.00001))),\
#                                     min(param[2][1]/4.0, param[2][3]+random.uniform(-0.0001, 0.0001)),  max(0.00001,param[2][4]+random.uniform(-0.0001, 0.0001))])
##)

#                if param[3][0] == "smooth":
#                     newPoint.append(["smooth", param[3][1]+random.uniform(-0.0001, 0.0001), max(0.0, min(1.1,param[3][2]+random.uniform(-0.0001, 0.00001)))])#

#                elif param[3][0] == "corner":
#                    newPoint.append(["corner", param[3][1]+random.uniform(-0.0001, 0.0001), max(0.0, min(1.1,param[3][2]+random.uniform(-0.0001, 0.00001))),\
#                                                          min(param[3][1]/4.0, param[3][3]+random.uniform(-0.0001, 0.0001)), max(0.1, param[3][4]+random.uniform(-0.0001, 0.0001))])
#                                                          #
#                individual.listOfDrawnPoints.append(newPoint)

#        individual.listOfDrawnPoints.append(reference[-1])

#        individual.level=None




def generateFunc(numgen, params):
    nbInd = params["nbInd"]
    pop = algorithm.Population()

    for i in range(nbInd):
        global reference
        pop.append(generateIndividual(reference))

    return pop

####################################################################################################
### Eval
###
def evaluationFunc(pop):

    print("Evaluation Function "+str(len(pop)))
    basedir=os.path.dirname(__file__)
    bestscore = -float(sys.maxint)

    filename=[]
    for ind in pop:

        fend =  "def evalField(x,y,z):\n   return shape.eval(primitives.Point(x,y,z))"
        f1 = crochetInFile.toPythonString(ind) + fend
        filename.append((f1, ind))

    #################### EXAMPLE USING THE SEQUENTIAL LAUNCHER #################################
    ### List of filename that contains the simulation to run
    scenefiles = ["scene.pyscn","controller.py", "shape.py"]
    filesandtemplates = []
    for scenefile in scenefiles:
        filesandtemplates.append( (open(basedir+"/"+scenefile).read(), scenefile) )

    runs = []
    for f1,ind in filename:
        runs.append( {"GENERATION": str(pop.id),
        "INDIVIDUAL": str(ind.id),
        "SHAPECONTENT": f1, "nbIterations":180,
        "LIBRARYPATH" : os.path.dirname(geneticalgo.__file__)
        } )
    results = launcher.startSofa(runs, filesandtemplates, launcher=launcher.SerialLauncher())

    for res in results:
        print("Results: ")
        print("    directory: "+res["directory"])
        print("        scene: "+res["scene"])
        print("      logfile: "+res["logfile"])
        print("     duration: "+str(res["duration"])+" sec")

    ### Associate the results to the individuals.
    for i in range(len(filename)):
        f1, ind = filename[i]
        ind.results = results[i]
        data = getJSONFragmentFrom( ind.results["logfile"] )


    if data == None:
        print "SOFA CRASHED DOWN!!  TRY TO LAUNCH MANUALLY scene.pyscn"
        ind.level = - float(sys.maxint)
    else:

        VerticalGap1 = data["VerticalGap1"]
        VerticalGap2 = data["VerticalGap2"]

        if VerticalGap1 < 0.0 or VerticalGap2 > 0.0:
            print "STRANGE BEHAVIOR, VerticalGap1 < 0.0 or VerticalGap2 > 0.0"
            level = VerticalGap1 + 10.0*VerticalGap2

        if abs(level) > 10.0 or VerticalGap1 == 0.0 or VerticalGap2 == 0.0:
            print "bad shape"
            ind.level = - float(sys.maxint)
        else:
            ind.level = level



    ## Wait here the results.
    for f1, ind in filename:
        if ind.level > bestscore:
            bestscore = ind.level

    return bestscore


####################################################################################################
### Select
###
def selectionFunc(pop, params):

    nbInd=params["nbInd"]

    if not isinstance(nbInd, int):
        raise ValueError, "nbInd is not given or not an int"
    pop.sortAccordingToLevel()
    pop.deleteFrom(params["nbInd"])
#    pop.sortAccordingToIndex()

    return pop
