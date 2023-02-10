import numpy as np
import re
import random
import logging

def random_position():
    RANDOM_RANGE_LOWER = -100
    RANDOM_RANGE_UPPER = 100
    return random.randint(RANDOM_RANGE_LOWER, RANDOM_RANGE_UPPER)

def random_boson_positions(nbosons):
    return [[random_position(), random_position(), random_position()] for _ in range(nbosons)]

def set_logger(logfile):
    logging.basicConfig(format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                        datefmt='%H:%M:%S',
                        level=logging.DEBUG,
                        handlers=[
                                            logging.FileHandler(logfile),
                                            logging.StreamHandler()
                                    ])

    global logger
    logger = logging.getLogger(__name__)
    return logger