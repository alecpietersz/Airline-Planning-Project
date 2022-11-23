from openpyxl import *

class AirportPair:
    def __init__(self,origin, destination, distance, demand):
        self.From       = origin
        self.To         = destination
        self.Distance   = distance
        self.Demand     = demand

class Airport:
    def __init__(self,name, icao, latitude, longitude, runway_length, population, gdp, hub):
        self.Name           = name
        self.ICAO           = icao
        self.Lat            = latitude
        self.Long           = longitude
        self.RunwayLength   = runway_length
        self.Population     = population
        self.GDP            = gdp        
        self.Hub            = hub



wb = load_workbook("Group_16_Aircraft_info.csv", read_only=True)
List_aircraft_info = tuple(wb["Group_16_Aircraft_info"].iter_rows())

wb = load_workbook("Group_16_Airport_info.csv", read_only=True)
List_our_airports = tuple(wb["Group_16_Airport_info"].iter_rows())

wb = load_workbook("Group_16_Demand.csv", read_only=True)
List_demand_pairs = tuple(wb["Group_16_Demand"].iter_rows())

wb = load_workbook("Group_16_Distances.csv", read_only=True)
List_demand_pairs = tuple(wb["Group_16_Distances"].iter_rows())

wb = load_workbook("Group_16_Annual_growth.csv", read_only=True)
List_demand_pairs = tuple(wb["Group_16_Annual_growth"].iter_rows())


for i, origin in enumerate(List_airport_demands[1:]):
    for j, destination in enumerate(origin[1:]):
        new_pair = AirportPair(List_airport_demands[i+1][0].value,List_airport_demands[0][j+1].value,int(List_airport_distances[i+1][j+1].value),int(destination.value))
        AirportPairs.append(new_pair)

for (name,hub) in List_airports[1:]:
    new = Airport(name.value,int(hub.value))
    Airports.append(new)

    