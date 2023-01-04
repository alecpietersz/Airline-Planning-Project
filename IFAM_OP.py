# Group 16
# Mauryze Brug (4700651)
# Alec Pietersz (5020328)
# Emma Zadeits (4671880)

from gurobipy import Model, quicksum, GRB, LinExpr
from numpy import *
from openpyxl import load_workbook
from time import *
import numpy as np
import itertools
import pickle
import xlsxwriter
import datetime
import InputGenOptions

def IFAM_Problem (L,P,N,K,RR,Arcs,Nodes,NGk, Delta, LF, LO, PO, solve):
    
    def ConstructModel(columns):

        model = Model("IFAM_CG") # LP model (this is an object)

        obj = LinExpr()

        
        # create decision variables
        t = {}
        for (itin1,itin2) in columns:            
            t[itin1.ID,itin2.ID] = model.addVar(vtype = 'I',name = f"t{itin1.ID}-{itin2.ID}")
            obj += t[itin1.ID,itin2.ID]*(itin1.Fare-RR[itin1.ID,itin2.ID]*itin2.Fare)
        
        f = {}
        for flight in L:
            for aircraft in K:
                f[flight.FN,aircraft.Type] = model.addVar(vtype = 'B', name =f"f{flight.FN}-{aircraft.Type}")
                obj += f[flight.FN,aircraft.Type]*flight.ACCost[aircraft.Type]

        y = {}
        for id, arc in Arcs.items():
            if arc['arc_type'] == 'ground':
                y[id,arc['ac_type']] = model.addVar(vtype='I',name = f"y{id}-{arc['ac_type']}")

        Z = {}
        for itin in PO:
            Z[itin.ID] = model.addVar(vtype='B',name = f"Z{itin.ID}")
            obj+= (1-Z[itin.ID])*itin.Demand*itin.Fare

        model.setObjective(obj,GRB.MINIMIZE)
        model.update()
        
        # generate constraints

        C1 = {}
        for flight in LF:
            C1[flight.FN] = model.addLConstr(quicksum(f[flight.FN,k.Type] for k in K),'=',1, name = f"C1{flight.FN}")

        C1B = {}
        for flight in LO:
            C1B[flight.FN] = model.addLConstr(quicksum(f[flight.FN,k.Type] for k in K),'<=',1, name = f"C1{flight.FN}")

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
            sum = 0
            for itin in columns:
                if (itin[1],itin[0]) in t:
                    sum+=Delta[flight.FN,itin[0]]*RR[itin[0],itin[1]]*t[itin[1],itin[0]]

            C4[flight.FN] = model.addLConstr(quicksum(aircraft.Seats*f[flight.FN,aircraft.Type] for aircraft in K) + 
                                            quicksum(Delta[flight.FN,itin[0].ID]*t[itin[0].ID,itin[1].ID] for itin in columns)-sum,
                                            '>=', quicksum(Delta[flight.FN,itin3.ID]*itin3.Demand for itin3 in P), name = f"C4{flight.FN}")

        C5 = {}
        for p in P:
            sum = 0
            for r in P:
                if (p.ID,r.ID) in columns:
                    sum+= t[p.ID,r.ID]
            C5[p.ID] = model.addLConstr(sum, '<=', p.Demand, name = f"C5{p.ID}")

        C6 = {}
        for p in PO:
            for leg in p.Legs:
                C6[p.ID,leg] = model.addLConstr(Z[p.ID]-quicksum(f[leg,k.Type] for k in K),"<=",0, name = f"C6{p.ID}-{leg}")

        C7 = {}
        for p in PO:
            C7[p.ID] = model.addLConstr(Z[p.ID]-quicksum(quicksum(f[leg,k.Type] for k in K) for leg in p.Legs),"=>",1-len(p.Legs),name=f"C7{p.ID}")
        
        model.update()
        model.write("IFAM_OP.lp")
        
        return model

    columns = []

    for itin1 in P:
        for itin2 in P:
            if itin1.ID == itin2.ID or itin2.ID != 0:
                continue
            else:
                columns.append((itin1,itin2))

    model = ConstructModel(columns)

    count = 0
    exit_condition = False

    while exit_condition is False and count<100:
        print("")
        print("START RELAXED ITERATION",count)
        exit_condition = True
        linear_relaxation = model.relax()
        # linear_relaxation.write(f"IFAM_OP_relaxed{count}.lp")
        linear_relaxation.optimize()

        Pi = {}
        Sigma = {}        
        for c in linear_relaxation.getConstrs():
            if c.ConstrName[1] == '4':
                Pi[c.ConstrName[2:]] = c.Pi
            elif c.ConstrName[1] == '5':
                Sigma[int(c.ConstrName[2:])] = c.Pi

        c = {}      
        for itin1 in P[1:]:
            for itin2 in P[1:]:
                c[itin1.ID,itin2.ID] = (itin1.Fare - quicksum(Pi[i] for i in itin1.Legs)) -RR[itin1.ID,itin2.ID]*(itin2.Fare-quicksum(Pi[j] for j in itin2.Legs))-Sigma[itin1.ID]
                
                if c[itin1.ID,itin2.ID].getValue() < 0:                    
                    if (itin1,itin2) not in columns:
                        columns.append((itin1,itin2))
                        print("COLIUMN ADDED",(itin1.ID,itin2.ID))
                        exit_condition = False
        
        model = ConstructModel(columns)

        count+=1

    print("")
    print("SOLVE NON-RELAXED")
    model = ConstructModel(columns)
    model.optimize()
    model.write("IFAM_OP_Final.sol") 
   
            
    print
    
    print ("------------------------------------------------------------------------")


if __name__ == '__main__':  

    start_time = time()

    file = open('input_data_options.pickle', 'rb')
    data = pickle.load(file)

    L,P,N,K,RR,Arcs,Nodes,NGk,Delta, LF, LO, PO = data

    # for i in P:
    #     print(vars(i))

    IFAM_Problem(L,P,N,K,RR,Arcs,Nodes,NGk, Delta, LF, LO, PO, solve=True)
    
    elapsed_time = time() - start_time

    print ("Run Time = ", elapsed_time)
    print ("END")