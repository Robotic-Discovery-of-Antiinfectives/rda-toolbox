#!/usr/bin/env python3

import numpy as np

def zfactor(positive_controls, negative_controls):
    return 1 - (3*(np.std(positive_controls) + np.std(negative_controls)) /
                abs(np.mean(positive_controls - np.mean(negative_controls))))

def minmax_normalization(x, minimum, maximum):
    return ((x - minimum) / (maximum - minimum)) * 100

def max_normalization(x, maximum):
    return (x / maximum) * 100
