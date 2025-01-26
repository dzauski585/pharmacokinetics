import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
import leanbodymass as lbm
import schnider_params as sp

weight = 75 #kg
height = 170 #cm
age = 85 #yr
gender = 'f'
model = 'schneider'



if model == 'schneider':
    params = sp.schnider_params_calc(age, weight, height)

print(params)

v1, v2, v3, v4, c5, v6, v7, v8 = params

print(v1)