from gurobipy import Model, quicksum, GRB
from numpy import *
from openpyxl import *
from time import *
import numpy as np
import itertools
import pickle
from RouteBased import FN_Problem

with open('routebased.pickle') as f:
    (AirportPairs, Aircrafts, Airports, fuel_price, block_time, energy_price, load_factor, valid_routes, routes2, routes1) = pickle.load(f)

FN_Problem(AirportPairs, Aircrafts, Airports, fuel_price, block_time, energy_price, load_factor, valid_routes, routes2, routes1, solve=False)