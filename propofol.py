import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
import leanbodymass as lbm

weight= 75 #kg
height= 170 #cm
age= 85 #yr
gender= 'f'
lbm= lbm.lbm_calc(gender, weight, height)


f = weight * 7 / 60 #administration rate weight* mg/Kg/h (/60 makes it /min)

c1_init = 0 #initial conc at central compartment (1)
c2_init = 0 #initial conc at rapid_peripheral compartment (2)
c3_init = 0 #initial conc at slow_peripheral compartment (3)

def schnider_params_calc(a,w,h,lbm):
    vd_central = 4.27 #L
    vd_rapid_peripheral = 18.9-0.391*(a-53) #L
    vd_slow_peripheral = 238 #L
    clearance_met = 1.89 + ((w - 77) * 0.0456) + ((lbm - 59) * -0.0681) + ((h - 177) * 0.0264)
    clearance_rapid_periph = 1.29-0.024 * (a - 53)
    clearance_slow_periph = 0.836 #all clearances L min-1
    
    #to make integration easier, I calculate elim constants
    k10 = clearance_met / vd_central
    k12 = clearance_rapid_periph / vd_central
    k21 = clearance_rapid_periph /vd_rapid_peripheral
    k13 = clearance_slow_periph / vd_central
    k31 = clearance_slow_periph / vd_slow_peripheral
    
    schneider_params=(vd_central,vd_rapid_peripheral,vd_slow_peripheral,k10,k12,k21,k13,k31)
    return schneider_params

vd_1,vd_2,vd_3,k10,k12,k21,k13,k31 = schnider_params_calc(age,weight,height,lbm)


timespan  = 1 * 60 #n hours, in minutes
t = np.linspace(0,timespan,100)
timesteps = timespan / 1000

c1= np.zeros(t.shape) #concentration at every t (mg/L), central compartment
c2= np.zeros(t.shape) #concentration at every t (mg/L), RapidPeriph compartment
c3= np.zeros(t.shape) #concentration at every t (mg/L), SlowPeriph compartment

#initial conditions, differential equations interactions
def deriv(y, t):
    c1, c2, c3 = y
    dc1dt = (f / vd_1) + c2 * k21 * vd_2 / vd_1 + c3 * k31 * vd_3 / vd_1 - c1 * (k10 + k12 + k13)
    dc2dt = c1 * k12 * vd_1  /vd_2 - c2 * k21
    dc3dt = c1 * k13  *vd_1 / vd_3 - c3 * k31
    
    return dc1dt, dc2dt, dc3dt

y0= c1_init, c2_init, c3_init

    #integrate pk equations over time grid t differential equaions
ret = odeint(deriv, y0, t)

c1_integrated, c2_integrated, c3_integrated = ret.T




fig= plt.figure('Schnider propofol', figsize=(14,10))
ax= fig.add_subplot(111)

l1, = ax.plot(t, c1_integrated,color='C1', linewidth= 2, label= 'Central')
l2, =ax.plot(t, c2_integrated,color='C2', ls=':',linewidth= 2, label= 'Rapid peripheral')
l3, =ax.plot(t, c3_integrated,color='C3', ls=':',linewidth= 1, label= 'Slow peripheral')

ax.set_ylim([0,4])
ax.legend(loc='best')
ax.set_xlabel('Time in minutes')
ax.set_ylabel("Concentration in micrograms/ml")
ax.set_title('Schenider propofol model plasma concentration')

plt.show()