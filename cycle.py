from typing import TYPE_CHECKING
import geopandas as gpd
from numpy import random
import osmnx as ox
import numpy as np
import chaospy

# make a random generator (0-1, uniform)
rng = np.random.default_rng()

def make_cycle(edges, elevations, idle_normal, a_kernel, dcc_kernel, speedkmh):
        
    ### randomizex stuff ###
    p_stop = 0.1 # stopping probability
    stops = 0 # how many stops on route? Start a counter.

    # I. Stop or not?
    edges_with_stops=[]
    edge_adder=0 # if no stop, add previous edge distance to edge_adder
    for edge in edges:
        # stop or not?
        if rng.random() <= p_stop:
            edges_with_stops.append(edge+edge_adder)
            edge_adder=0
            stops = stops + 1
        else:
            edge_adder = edge_adder + edge
    
    # If there are no stops, make the route into one edge.
    if not edges_with_stops:
        edges_with_stops.append(sum(edges))

    t_data=[0]
    v_data=[0]

    v_max = speedkmh/3.6 # vehicle speed m/s
    ts=1 # timestep
    a=0.5 # m/s^2

    v=0 # speed is zero in the beginning
    d=0 # distance is zero in the beinning
    t=0 # time is zero in the beginning
    d_tot=0
    t_i=0

    a_data = []
    dcc_data= []

    # 2. one edge at a time
    for edge in edges_with_stops:
        dcc = -dcc_kernel.sample(1)[0]
        t_idle = idle_normal.sample(1)[0]
        a = a_kernel.sample(1)[0]
        dcc_data.append(dcc)
        a_data.append(a)

        # acceleration (third of the edge)
        while d < edge/3:
            # Increase speed if not max
            if v < v_max:
                v = v + a*ts
            d = d + v*ts
            t=t+ts
            t_data.append(t)
            v_data.append(v)

        # steady state and deceleration
        # 1. calculate time for deceleration
        t_stop = v/dcc
        d_stop = v*t_stop - 0.5*dcc*t_stop**2
        # 2. cruise while you can
        while d < (edge-d_stop):
            d = d + v*ts
            t=t+ts
            t_data.append(t)
            v_data.append(v)
        # 3. decelerate    
        while t_i < t_stop:
            if v > 0:
                v = v - dcc*ts
                if v < 0: # don't let v go below zero
                    v=0
            d = d + v*ts
            t_i=t_i+ts
            t=t+ts
            t_data.append(t)
            v_data.append(v)
        # 4. idle, possibly
        t_id=0
        while t_id < t_idle:
            t_id = t_id+ts
            t=t+ts
            t_data.append(t)
            v_data.append(0) #standing still
        t_i=0
        d_tot=d_tot+d
        d=0

    ### Defining road slope
    # There is always one more node (elevation point) than edges
    g_data = []
    for i in range(1,len(elevations)):
        dy = elevations[i]-elevations[i-1]
        dx = edges[i-1]
        theta = np.arctan(dy/dx) # road slope in degrees
        g_data.append(theta)
        #grade=10*((2*np.pi)/360) # 10 degs

    ele_data = [0]
    ele = 0
    d_e = 0
    d_e_tot = 0
    j=0 # edge iterations
    for i in range(0,len(v_data)-1):
        delta_de=v_data[i]*ts
        d_e = d_e + delta_de # elevation is dependent on distance
        d_e_tot= d_e_tot + delta_de
        if len(g_data) > j:
            ele_data.append(np.tan(g_data[j])*delta_de + ele)
            ele = ele + np.tan(g_data[j])*delta_de
        else:
            ele_data.append(ele)
        if len(edges) > j:
            if d_e > edges[j]:
                if len(edges) >= j:
                    d_e = 0 # zero the d_e for next edge
                    j=j+1

    return t_data, v_data, ele_data, sum(edges)-d_tot, stops, np.mean(a_data), np.mean(dcc_data)


def make_noise(v_data, sigma=0.1, ma_window=4):
    #SynthetiNoising the signal
    autoErr = chaospy.Normal(mu=0, sigma=sigma)
    autoCorr = 0.95
    v_data_ar =[]
    for i in range(1,len(v_data)+1):
        v_data_ar.append(autoCorr*v_data[i-1] + autoErr.sample(1)[0])

    moving_averages = []
    window_size = ma_window
    for i in range(0,len(v_data_ar)):
        if len(v_data_ar) - window_size <= 0:
            window_size -= 1
            print(window_size)
        window = v_data_ar[i : i + window_size]
        window_average = sum(window) / window_size
        moving_averages.append(window_average)
        i += 1

    if ma_window > 3:
        # Finalizing
        #Beginning and end-ramps, never go below 0
        moving_averages[0] = 0
        moving_averages[-1] = 0
        window_size = 2
        for i in range(0, len(moving_averages)):
            if i < 15:
                window = moving_averages[i : i + window_size]
                moving_averages[i] = sum(window) / window_size
            if len(moving_averages)-i < 15:
                if len(v_data_ar) - window_size <= 0:
                    window_size -= 1
                window = moving_averages[i : i + window_size]
                moving_averages[i] = sum(window) / window_size

            if moving_averages[i] < 0:
                moving_averages[i] = 0

        moving_averages[0] = 0
        moving_averages[-1] = 0
    else:
        moving_averages[-2] = moving_averages[-3]
        moving_averages[-1] = moving_averages[-2]

    return moving_averages