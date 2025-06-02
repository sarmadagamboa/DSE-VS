import numpy as np

def calc_payload_pointing(baffle_diameter=0.025, mirror_distance=2):

    req_angle = np.arctan((baffle_diameter/2)/mirror_distance)

    return req_angle

print(calc_payload_pointing())