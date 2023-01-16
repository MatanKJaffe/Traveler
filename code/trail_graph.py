import osmnx
from pyrosm import OSM, get_data
import geopandas as gpd
import pandas as pd
from sklearn.neighbors import BallTree
import numpy as np
import mapclassify as mc
import matplotlib.pyplot as plt
import time
import networkx as nx
import igraph as ig
import osmnx as ox
import folium
import json




osm = OSM(get_data("Southern California"))

graph_type = {}

n_drive,  e_drive = osm.get_network(nodes=True, network_type="driving")
n_cycling, e_cycling   = osm.get_network(nodes=True, network_type="cycling")
n_walk,   e_walk       = osm.get_network(nodes=True, network_type="walking")
n_service, e_service   = osm.get_network(nodes=True, network_type="driving+service")

graph_type['drive']     = ox.add_edge_travel_times(ox.add_edge_speeds(osm.to_graph(n_drive, e_drive,  graph_type="networkx")), extra_kwargs={"hv":{"car":120}})
graph_type['walk']   = ox.add_edge_travel_times(ox.add_edge_speeds(osm.to_graph(n_walk, e_walk, graph_type="networkx")))
graph_type['cycle']   = ox.add_edge_travel_times(ox.add_edge_speeds(osm.to_graph(n_cycling, e_cycling, graph_type="networkx")))
graph_type['service']   = ox.add_edge_travel_times(ox.add_edge_speeds(osm.to_graph(n_service, e_service, graph_type="networkx")))



def get_route(source_geo, dest_geo, go_type='drive', weight='travel_time',plot=True):
    source_node = ox.get_nearest_node(graph_type[go_type], source_geo)
    target_node = ox.get_nearest_node(graph_type[go_type], dest_geo)

    route = nx.shortest_path(graph_type[go_type], source_node, target_node, weight=weight)
    
    edge_lengths = ox.utils_graph.get_route_edge_attributes(graph_type[go_type], route, 'length') 
    edge_travel_time = ox.utils_graph.get_route_edge_attributes( graph_type[go_type], route, 'travel_time') 
    total_route_length = round(sum(edge_lengths), 1)
    route_travel_time  = round(sum(edge_travel_time)/60, 2)
    if plot:
        ox.plot_graph_route(graph_type[go_type], route, node_size=0, figsize=(40,40))
    return route, total_route_length, route_travel_time

def main():
    hospitals = json.loads("hospitals.json") # read from api hospitals list
    source_loc = (34.018255, -118.313290) # just for test you must get it from gps
    nearest_hospitals = []
    for hospital in hospitals:
        hospital_loc = (hospital['lat'], hospital['long'])
        hospital_score = hospital['score'] # this is score of hospital in api just for rating
        route, route_lenght, route_time = get_route(source_loc, hospital_loc)
        # route_lenght -> distance between us and hospital in meter
        # route_time   -> duration of startpoint to hospital in minutes
        # route        -> array of point by point to hospital 
        hospital_score -= (route_time * 5) # dec 5 score per minute to arrive at the hospital
        nearest_hospitals.append({'id':hospital['id'],'score': hospital_score})
    newlist = sorted(nearest_hospitals, key=lambda d: d['score'], reverse=True) 
    print(newlist)
    # new list is most populer and nearst hospital from start point
    from markupsafe import Markup, escape