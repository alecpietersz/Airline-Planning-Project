# Group 16
# Mauryze Brug (4700651)
# Alec Pietersz (5020328)
# Emma Zadeits (4671880)

from gurobipy import Model, quicksum, GRB
from numpy import *
from openpyxl import load_workbook
from time import *
import numpy as np
import itertools
import pickle
import xlsxwriter
import datetime

def IFAM_Problem (AirportPairs, Aircrafts, Airports, fuel_price, block_time, energy_price, load_factor, Routes, Routes2, Routes1, solve):
    
    model = Model("IFAM") # LP model (this is an object)

    # create decision variables    
   
    model.update()
    
    # generate constraint

   

    # run model
    if solve:
        model.update()
        model.write("IFAM.lp")
        model.optimize()
        model.write("IFAM.sol")
    else:
        model.update()
        model.read("IFAM.sol")
        model.setParam('TimeLimit', 10)
        model.optimize()

    status = model.status
    if status != GRB.Status.OPTIMAL:
        if status == GRB.Status.UNBOUNDED:
            print('The model cannot be solved because it is unbounded')
        elif status == GRB.Status.INFEASIBLE:
            print('The model is infeasible; computing IIS')
            model.computeIIS()
            print('\nThe following constraint(s) cannot be satisfied:')
            for c in model.getConstrs():
                if c.IISConstr:
                    print('%s' % c.constrName)
        elif status != GRB.Status.INF_OR_UNBD:
            print('Optimization was stopped with status %d' % status)

    print

            
    print
    print ("Objective Function =", model.ObjVal/1.0)
    print ("------------------------------------------------------------------------")


if __name__ == '__main__':  

    start_time = time()

    file = open('input_data.pickle', 'rb')
    data = pickle.load(file)

    L,P,N,K,RR,Arcs,Nodes = data

    IFAM_Problem(L,P,N,K,RR,Arcs,Nodes, solve=True)
    
    elapsed_time = time() - start_time

    print ("Run Time = ", elapsed_time)
    print ("END")