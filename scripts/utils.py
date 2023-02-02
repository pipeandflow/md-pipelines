import numpy as np
import re
import random

def random_position():
    RANDOM_RANGE_LOWER = -100
    RANDOM_RANGE_UPPER = 100
    return random.randint(RANDOM_RANGE_LOWER, RANDOM_RANGE_UPPER)

def random_boson_positions(nbosons):
    return [[random_position(), random_position(), random_position()] for _ in range(nbosons)]
