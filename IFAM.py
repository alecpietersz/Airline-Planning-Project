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
import InputGen

def IFAM_Problem (L,P,N,K,RR,Arcs,Nodes,NGk, Delta, solve):
    
    model = Model("IFAM") # LP model (this is an object)

    # create decision variables
    t = {}
    for itin1 in P:
        for itin2 in P:
            t[itin1.ID,itin2.ID] = model.addVar(obj=(itin1.Fare-RR[itin1.ID,itin2.ID]*itin2.Fare), vtype = 'I',name = f"t{itin1.ID}-{itin2.ID}")
    
    f = {}
    for flight in L:
        for aircraft in K:
            f[flight.FN,aircraft.Type] = model.addVar(obj=flight.ACCost[aircraft.Type], vtype = 'B', name =f"f{flight.FN}-{aircraft.Type}")

    y = {}
    for id, arc in Arcs.items():
        if arc['arc_type'] == 'ground':
            y[id,arc['ac_type']] = model.addVar(obj=0,vtype='I',name = f"y{id}-{arc['ac_type']}")    

    model.update()
    
    # generate constraints

    C1 = {}
    for flight in L:
        C1[flight.FN] = model.addLConstr(quicksum(f[flight.FN,k.Type] for k in K),'=',1, name = f"C1{flight.FN}")

    C2 = {}
    for key, node in Nodes.items():
        C2[key] = model.addLConstr(y[node['g_o_arc'],key[0]]+ quicksum(f[Arcs[i]['arc_type'],key[0]] for i in node['o_arcs'])-
                                  y[node['g_e_arc'],key[0]]- quicksum(f[Arcs[i]['arc_type'],key[0]] for i in node['e_arcs']),
                                  '=', 0, name = f"C2{key}".replace(" ", ""))
    
    C3 = {}
    for aircraft in K:
        C3[aircraft.Type] = model.addLConstr(quicksum(y[a,aircraft.Type] for a in NGk[aircraft.Type] if Arcs[a]['arc_type'] == 'ground') + quicksum(f[Arcs[a]['arc_type'],aircraft.Type] for a in NGk[aircraft.Type] if Arcs[a]['arc_type'] != 'ground'),'<=',aircraft.Units, name = f"C3{aircraft.Type}")

    C4 = {}
    for flight in L:
        C4[flight.FN] = model.addLConstr(quicksum(aircraft.Seats*f[flight.FN,aircraft.Type] for aircraft in K) + 
                                        quicksum(quicksum(Delta[flight.FN,itin1.ID]*t[itin1.ID,itin2.ID] for itin1 in P)for itin2 in P)-
                                        quicksum(quicksum(Delta[flight.FN,itin1.ID]*RR[itin1.ID,itin2.ID]*t[itin1.ID,itin2.ID] for itin1 in P)for itin2 in P),
                                        '>=', quicksum(Delta[flight.FN,itin3.ID]*itin3.Demand for itin3 in P), name = f"C4{flight.FN}")

    C5 = {}
    for p in P:
        C5[p.ID] = model.addLConstr(quicksum(t[p.ID,r.ID]for r in P), '<=', p.Demand, name = f"C5{p.ID}")

    C6 = {}
    for p in P:
        C6[p.ID] = model.addLConstr(t[p.ID,p.ID], '=', 0, name = f"C6{p.ID}")

    # run model
    
    model.update()
    model.write("IFAM.lp")
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
    else:
        model.write("IFAM.sol")
        print ("Objective Function =", model.ObjVal/1.0)

    print

    
            
    print
    
    print ("------------------------------------------------------------------------")


if __name__ == '__main__':  

    start_time = time()

    file = open('input_data.pickle', 'rb')
    data = pickle.load(file)

    L,P,N,K,RR,Arcs,Nodes,NGk,Delta = data

    IFAM_Problem(L,P,N,K,RR,Arcs,Nodes,NGk, Delta, solve=True)
    
    elapsed_time = time() - start_time

    print ("Run Time = ", elapsed_time)
    print ("END")