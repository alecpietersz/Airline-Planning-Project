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

class Flight:
    def __init__(self, fn, origin, destination, dep_time, arr_time, ac_cost):
        self.FN             = fn
        self.Origin         = origin
        self.Destination    = destination
        self.DepTime        = dep_time
        self.ArrTime        = arr_time
        self.ACCost         = ac_cost

class Itinerary:
    def __init__(self, id, origin, destination, demand, fare, stops, legs):
        self.ID             = id
        self.Origin         = origin
        self.Destination    = destination
        self.Demand         = demand
        self.Fare           = fare
        self.Stops          = stops
        self.Legs           = legs

class Aircraft:
    def __init__(self, type, units, seats, tat):
        self.Type  = type
        self.Units = units
        self.Seats = seats
        self.TAT   = tat

L      = []
P      = []
N      = []
K      = []
RR     = {}

wb = load_workbook("Group_16_copy.xlsx", read_only=True)
List_flights            = tuple(wb["Flight"].iter_rows())
List_itineraries        = tuple(wb["Itinerary"].iter_rows())
List_recapture_rate     = tuple(wb["Recapture Rate"].iter_rows())
List_aircraft           = tuple(wb["Aircraft"].iter_rows())

for flight in List_flights[1:]:
    if flight[0].value:
        ac_cost = {}
        ac_cost[List_flights[0][5].value] = flight[5].value
        ac_cost[List_flights[0][6].value] = flight[6].value
        ac_cost[List_flights[0][7].value] = flight[7].value
        ac_cost[List_flights[0][8].value] = flight[8].value
        ac_cost[List_flights[0][9].value] = flight[9].value
        new_flight = Flight(flight[0].value,flight[1].value,flight[2].value,flight[3].value,flight[4].value,ac_cost)
        L.append(new_flight)

for itinerary in List_itineraries[1:]:
    if itinerary[0].value:
        legs = []
        legs.append(itinerary[6].value)
        if itinerary[7].value != 0:
            legs.append(itinerary[7].value)
        new_itinerary = Itinerary(itinerary[0].value,itinerary[1].value,itinerary[2].value,itinerary[3].value,itinerary[4].value,itinerary[5].value,legs)
        P.append(new_itinerary)

for rr_pair in List_recapture_rate[1:]:
    RR[rr_pair[0].value,rr_pair[1].value] = rr_pair[2].value

for aircraft in List_aircraft[1:]:
    if aircraft[0].value:
        new_aircraft = Aircraft(aircraft[0].value,aircraft[1].value,aircraft[2].value,aircraft[3].value)
        K.append(new_aircraft)

for flight in K:
    print(vars(flight))

# print(RR.values())
