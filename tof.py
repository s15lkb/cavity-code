# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 12:01:17 2017

@author: valentin.metillon
"""



positions = {'BAC': 0, 'R1' : 56.8, 'C1': 101.9, 'R2':147, 'C2': 192.1,'R3': 237.2, 'D1': 278, 'D2': 308}

def tof(pos1, pos2, v=250):
    #returns the time of flight between position 1 and position 2 in microseconds. 
    #pos1 and pos2 should be positions expressed in millimeters or keys to the positions dictionary
    if type(pos1)==str:
        pos1 = positions[pos1]
    if type(pos2)==str:
        pos2 = positions[pos2]
    return((pos2-pos1)*1000/v)

