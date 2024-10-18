import sys

def sim_volume(bug_N, bug_M, bug_spacing, bug_diameter, box_height):
    box_len = bug_N*(bug_spacing*bug_diameter)
    box_wid = bug_M*(bug_spacing*bug_diameter)
    box_vol = box_len*box_wid*box_height
    return(box_vol)
