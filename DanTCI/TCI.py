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


leanbm = lbm.lbm_calc(gender, weight, height)

if model == 'schneider':
    params = sp.schnider_params_calc(age, weight, height, leanbm)

print(params)