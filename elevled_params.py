import numpy as np
import fatfreemass as ffm
'''
D.J. Eleveld, P. Colin, A.R. Absalom, M.M.R.F. Struys
Pharmacokinetic–pharmacodynamic model for propofol for broad application in anaesthesia and sedation
Br J Anaesth, 120 (2018), pp. 942-959

https://www.bjanaesthesia.org/article/S0007-0912(18)30051-5/fulltext
'''



#%% Model parameters
t1 = 6.28        #L
t2 = 25.5        #L
t3 = 273          #L
t4 = 1.79        #L/min
t5 = 1.75        #L/min
t6 = 1.11        #L/min
t7 = 0.191
t8 = 42.30       #weeks
t9 = 9.06
t10 = -0.0156
t11 = -0.00286
t12 = 33.6       #Kg
t13 = -0.0138
t14 = 68.3       #weeks
t15 = 2.10       #L/min
t16 = 1.30 
t17 = 1.42 
t18 = 0.68
e1 = 0.610       #v1
e2 = 0.565       #v2
e3 = 0.597       #v3
e4 = 0.265       #Cl
e5 = 0.346       #Q2
e6 = 0.209       #Q3
e7 = 0.463       #Residual error
#%% without etas

#%%
def eleveld_params_calc(g, a, w, h):
    
    f_ageing = lambda x: np.exp(x * (a - 35))
    f_sigmoid = lambda x, E50, l: x**l/(x**l + E50**l)

    f_central = lambda x: f_sigmoid(x, t12, 1)
    f_cl_maturation =  f_sigmoid(a * 52 + 40 / 52, t8, t9)
    f_cl_maturation_ref = f_sigmoid(35 * 52 + 40, t8, t9)
    quot_cl_mat = f_cl_maturation / f_cl_maturation_ref
    
    f_q3_maturation =  f_sigmoid(a * 52 + 40, t14, 1)
    f_q3_maturation_ref = f_sigmoid(35 * 52 + 40, t14, 1)
    quot_f_q3_mat = f_q3_maturation / f_q3_maturation_ref
    
    f_opiates = lambda x: np.exp(x * a) #assume always opioids, otherwise 1
    
    v1  = t1 * (f_central(w) / f_central(70)) #vd_central
    v2  = t2 * f_ageing(t10) * (w / 70) #vd_rapid_peripheral
    v3  = t3 * ffm.fatfreemass(g, a, w, h) / ffm.fatfreemass('m', 35, 70, 170) * f_opiates(t13) #vd_slow_peripheral
    
    if g == 'm':
        theta = t4
    elif g == 'f':
        theta = t15
    cl = theta * ((w / 70)**0.75) * quot_cl_mat * f_opiates(t11)
    q2 = t5 * ((v2 / 25.5)**0.75) * (1 + t6 * (1 - f_q3_maturation))
    q3 = t6 * ((v3 / 273)**0.75) * quot_f_q3_mat 
    
    k10 = cl / v1
    k12 = q2 / v1
    k21 = q2 / v2
    k13 = q3 / v1
    k31 = q3 / v3
    
    Ce50 = 3.08 * f_ageing(-0.00635)
    ke0 = 0.146 * (w / 70)**(-0.25)
    
    model= 'Eleveld'
    drug_name= 'Propofol'
    units= (r'$\mu g\ ml^{-1}$', r'mg$\ Kg^{-1}\ h^{-1}$')
    
    
    params=[drug_name, 
            model,
            units, Ce50,
            v1,
            v2,
            v3,
            k10,
            k12,k21,
            k13,k31,
            ke0
            ]
    return params