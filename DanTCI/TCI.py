import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
import leanbodymass as lbm

weight= 75 #kg
height= 170 #cm
age= 85 #yr
gender= 'f'

lbm = lbm.lbm_calc(gender, weight, height)
