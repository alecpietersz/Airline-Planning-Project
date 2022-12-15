
from gurobipy import Model, quicksum, GRB
from numpy import *
from openpyxl import load_workbook
from time import *
import numpy as np
import itertools
import pickle

class AirportPair:
    def __init__(self,origin, destination, distance, demand2030):
        self.From       = origin
        self.To         = destination
        self.Distance   = distance
        self.Demand2030 = demand2030

class DemandAirportPair:
    def __init__(self,origin, destination, demand2020):
        self.From       = origin
        self.To         = destination
        self.Demand2020 = demand2020

class Airport:
    def __init__(self,name, icao, latitude, longitude, runway_length, population2020, population2030, gdp, hub):
        self.Name           = name
        self.ICAO           = icao
        self.Lat            = latitude
        self.Long           = longitude
        self.Runway         = runway_length
        self.Population2020 = population2020
        self.Population2030 = population2030
        self.GDP            = gdp        
        self.Hub            = hub

class Aircraft:
    def __init__(self,name, speed, seats, tat, charging, ac_range, runway, lease_cost, op_cost, time_cost, fuel_cost, energy_cost):
        self.Name           = name
        self.Speed          = speed
        self.Seats          = seats
        self.TAT            = tat/60
        self.Charging       = charging/60
        self.Range          = ac_range
        self.Runway         = runway
        self.Lease          = lease_cost
        self.OperatingCost  = op_cost
        self.TimeCost       = time_cost
        self.FuelCost       = fuel_cost
        self.EnergyCost     = energy_cost
    
class Route:
    def __init__(self, airports, distance, delta, precedent, subsequent, id, legs):
        self.ID         = id
        self.Airports   = airports
        self.Distance   = distance
        self.Delta      = delta
        self.Precedent  = precedent
        self.Subsequent = subsequent        
        self.Legs       = legs

# AirportPairs, Aircrafts, Airports, fuel_price, block_time, energy_price, load_factor, valid_routes, routes2, routes1
def FN_Problem (AirportPairs, Aircrafts, Airports, fuel_price, block_time, energy_price, load_factor, Routes, Routes2, Routes1, solve):
    
    model = Model("FN")                # LP model (this is an object)
    
    w = {}
    for r in range(len(Routes)):
        for n in range(len(Routes)):                              # Decision Variables (DVs
            for m in range(len(AirportPairs)):
                AirPair = AirportPairs[m]
                if AirPair.Distance > 0:
                    yield_rpk = 5.9*AirPair.Distance**-0.76+0.043
                else:
                    yield_rpk = 0     
                w[AirPair.From,AirPair.To,r,n] = model.addVar(obj=-yield_rpk*AirPair.Distance*0.9, vtype ="I",
                                                        name = ("w"+AirPair.From+'-'+AirPair.To+'-r'+str(Routes[r].ID)+'-n'+str(Routes[n].ID)).replace(" ", ""))

    x = {}
    for r in range(len(Routes)):                              # Decision Variables (DVs)
        for m in range(len(AirportPairs)):
            AirPair = AirportPairs[m]
            if AirPair.Distance > 0:
                yield_rpk = 5.9*AirPair.Distance**-0.76+0.043
            else:
                yield_rpk = 0 
            x[AirPair.From,AirPair.To,r] = model.addVar(obj=-yield_rpk*AirPair.Distance, vtype ="I",
                                                    name = ("x"+AirPair.From+'-'+AirPair.To+'-r'+str(Routes[r].ID)).replace(" ", ""))

    z = {}
    for r in range(len(Routes)):                              # Decision Variables (DVs)
        for k in range(len(Aircrafts)):              
            acft = Aircrafts[k]
            cost = 0
            for l in Routes[r].Legs:
                AirPair = next(airp for airp in AirportPairs if airp.From == l[0] and airp.To == l[1])
                gi = next(airp for airp in Airports if airp.ICAO == l[0]).Hub
                gj = next(airp for airp in Airports if airp.ICAO == l[1]).Hub
                cost += (acft.OperatingCost + acft.TimeCost*AirPair.Distance/acft.Speed + acft.FuelCost*fuel_price/1.5*AirPair.Distance)*(1-(2- gi - gj)*0.3) +energy_price*acft.EnergyCost*AirPair.Distance/acft.Range
            z[r,k] = model.addVar(obj=cost, vtype ="I",
                                                    name = ("z-r"+str(Routes[r].ID)+acft.Name).replace(" ", ""))

    ACk = {}                              # Decision Variables (DVs)
    for k in range(len(Aircrafts)):
        acft = Aircrafts[k]
        ACk[k] = model.addVar(obj=acft.Lease, vtype ="I", name = ("AC"+acft.Name).replace(" ", ""))

    model.update()

    Demand1 = {}                     # build 'capacity' constraints
    for m in range(len(AirportPairs)):
        AirPair = AirportPairs[m]
        Demand1[AirPair.From,AirPair.To] = model.addLConstr(quicksum(x[AirPair.From,AirPair.To,r]+quicksum(w[AirPair.From,AirPair.To,r,n] for n in range(len(Routes))) for r in range(len(Routes))),
                                                             '<=', AirPair.Demand2030, name = 'dem1'+AirPair.From+'-'+AirPair.To)


    # [item for item in accounts if item.get('id')==2]
    # [airp for airp in Airports if airp.get('Name')==AirportPairs[m].From].Name
    # INDX is the airport for which name == AirportPairs[m].From

    Demand2 = {}                       # build 'capacity' constraints
    for m in range(len(AirportPairs)):
        for r in range(len(Routes)):           
            AirPair = AirportPairs[m]
            Demand2[AirPair.From,AirPair.To,r] = model.addLConstr(x[AirPair.From,AirPair.To,r],
                                                                '<=', AirPair.Demand2030 * Routes[r].Delta[AirPair.From,AirPair.To], 
                                                                name = 'dem2'+AirPair.From+'-'+AirPair.To+'-r'+str(Routes[r].ID))


    Demand3 = {}                       # build 'capacity' constraints
    for m in range(len(AirportPairs)):
        for r in range(len(Routes)):
            for n in range(len(Routes)):
                AirPair = AirportPairs[m]
                Demand3[AirPair.From,AirPair.To] = model.addLConstr(w[AirPair.From,AirPair.To,r,n],
                                                                    '<=', AirPair.Demand2030*Routes[r].Delta[AirPair.From,AirPair.To]*Routes[n].Delta[AirPair.From,AirPair.To], 
                                                                    name = 'dem3'+AirPair.From+'-'+AirPair.To+'r'+str(Routes[r].ID)+'n'+str(Routes[n].ID))

    H = next(airp for airp in Airports if airp.Hub == 0).ICAO

    Flow1 = {}                       # build 'capacity' constraints
    for r in range(len(Routes)):
        Flow1[r] = model.addLConstr(quicksum(x[H,m,r] for m in Routes[r].Subsequent[H])
                                    + quicksum(quicksum(quicksum(w[p.ICAO,m,n,r] for m in Routes[r].Subsequent[H]) for p in Airports) for n in range(len(Routes))),
                                    '<=', quicksum(z[r,k]*Aircrafts[k].Seats*load_factor for k in range(len(Aircrafts))), name = 'flow1r'+str(Routes[r].ID))

    Flow2 = {}
                           # build 'capacity' constraints
    for r in range(len(Routes)):
        if Routes[r] in Routes2:
            AirPair = AirportPairs[m]
            i = Routes[r].Subsequent[H][0]
            j = Routes[r].Subsequent[H][1]
            Flow2[r] = model.addLConstr(quicksum(x[i,m,r] for m in Routes[r].Subsequent[j])+quicksum(x[m,j,r] for m in Routes[r].Precedent[i])
                                        + quicksum(quicksum(w[p.ICAO,j,n,r] for n in range(len(Routes))) for p in Airports)+ 
                                        quicksum(quicksum(w[i,p.ICAO,r,n] for n in range(len(Routes))) for p in Airports),
                                        '<=', quicksum(z[r,k]*Aircrafts[k].Seats*load_factor for k in range(len(Aircrafts))), name = 'flow2r'+str(Routes[r].ID))
        # print('')
        # print(Routes[r].Airports)
        # print("i",i)
        # print("j",j)
        # print("Sj",Routes[r].Subsequent[j])
        # print("Pi",Routes[r].Precedent[i])
        


    Flow31 = {}                       # build 'capacity' constraints
    for r in range(len(Routes)):  
        if Routes[r] in Routes2: 
            i = Routes[r].Subsequent[H][1]        
            Flow31[r] = model.addLConstr(quicksum(x[m,H,r] for m in Routes[r].Precedent[i])
                                        + quicksum(quicksum(quicksum(w[m,p.ICAO,r,n] for m in Routes[r].Precedent[i]) for p in Airports) for n in range(len(Routes))),
                                        '<=', quicksum(z[r,k]*Aircrafts[k].Seats*load_factor for k in range(len(Aircrafts))), name = 'flow3r'+str(Routes[r].ID))

    Flow32 = {}                       # build 'capacity' constraints
    for r in range(len(Routes)): 
        if Routes[r] in Routes1:  
            i = Routes[r].Subsequent[H][0]        
            Flow32[r] = model.addLConstr(quicksum(x[m,H,r] for m in Routes[r].Precedent[i])
                                        + quicksum(quicksum(quicksum(w[m,p.ICAO,r,n] for m in Routes[r].Precedent[i]) for p in Airports) for n in range(len(Routes))),
                                        '<=', quicksum(z[r,k]*Aircrafts[k].Seats*load_factor for k in range(len(Aircrafts))), name = 'flow3r'+str(Routes[r].ID))



    #AirportPairs[m].Distance/Aircrafts[0].Speed+Aircrafts[0].LTO)*z[AirportPairs[m].From,AirportPairs[m].To]

    ACProductivity = {}                       # build 'capacity' constraints
    for k in range(len(Aircrafts)):
        acft = Aircrafts[k]
        tat_list = {}
        for r in range(len(Routes)):
            tat = 0
            for l in Routes[r].Legs:
                AirPair = next(airp for airp in AirportPairs if airp.From == l[0] and airp.To == l[1])
                gj = next(airp for airp in Airports if airp.ICAO == l[1]).Hub
                tat += acft.TAT*(1+0.5*(1-gj))
            tat_list[r] = tat
        # print(tat_list)
             
        ACProductivity[k] = model.addLConstr(quicksum((Routes[r].Distance/acft.Speed + acft.Charging + tat_list[r])*z[r,k] for r in range(len(Routes))) ,'<=', block_time*ACk[k], name = 'ACProductivity'+Aircrafts[k].Name)

    Runway = {}
    for k in range(len(Aircrafts)):
        for r in range(len(Routes)):
            acft = Aircrafts[k]
            rwy_list = []
            for apt in Routes[r].Airports:
                rwy = next(airp for airp in Airports if airp.ICAO == apt).Runway
                rwy_list.append(rwy)
            min_rwy = min(rwy_list)

            if min_rwy >= acft.Runway:
                a = 10000
            else:
                a = 0
            Runway[r,k] = model.addLConstr(z[r,k],'<=', a, name = 'runway'+AirPair.From+'-'+AirPair.To+acft.Name)

    Range = {}
    for k in range(len(Aircrafts)):
        for r in range(len(Routes)):
            
            if Routes[r].Distance <= Aircrafts[k].Range:
                a = 10000
            else:
                a = 0
            Range[r,k] = model.addLConstr(z[r,k],'<=', a, name = 'range'+AirPair.From+'-'+AirPair.To+Aircrafts[k].Name)

    if solve:
        model.update()
        # model.setParam('TimeLimit', 1*60)
        model.write("RB_Model.lp")
        model.optimize()
        model.write("RB.sol")
    else:
        model.update()
        model.read("RB.sol")
        model.setParam('TimeLimit', 1)
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
        # exit(0)

    print
    
    
    for r in range(len(Routes)):
        for m in range(len(AirportPairs)):
            if  x[AirportPairs[m].From,AirportPairs[m].To,r].X>0:
                    print(x[AirportPairs[m].From,AirportPairs[m].To,r].VarName,x[AirportPairs[m].From,AirportPairs[m].To,r].X)
                    for n in range(len(Routes)):
                        if  w[AirportPairs[m].From,AirportPairs[m].To,r,n].X>0:
                         print(w[AirportPairs[m].From,AirportPairs[m].To,r,n].VarName,w[AirportPairs[m].From,AirportPairs[m].To,r,n].X)
                    
        for k in range(len(Aircrafts)):
            if z[r,k].X >0:
                print(z[r,k].VarName,z[r,k].X)
    
    for k in range(len(Aircrafts)):
        if ACk[k].X > 0:
            print(ACk[k].VarName,ACk[k].X)

            
    print
    print ("Objective Function =", model.ObjVal/1.0)
    print ("------------------------------------------------------------------------")


if __name__ == '__main__':
    

    start_time = time()
    # RUN MCF PROBLEM

    file = open('routebased.pickle', 'rb')
    data = pickle.load(file)

    AirportPairs, Aircrafts, Airports, fuel_price, block_time, energy_price, load_factor, valid_routes, routes2, routes1 = data

    FN_Problem(AirportPairs, Aircrafts, Airports, fuel_price, block_time, energy_price, load_factor, valid_routes, routes2, routes1, solve=True)
    
    elapsed_time = time() - start_time

    print ("Run Time = ", elapsed_time)
    print ("END")