'''
Schnider, Thomas W. Dr med; Minto, Charles F. MB, ChB; Gambus, Pedro L. MD; Andresen, Corina MD; 
Goodale, David B. DDS, PhD; Shafer, Steven L. MD; Youngs, Elizabeth J. MD. 
The Influence of Method of Administration and Covariates on the Pharmacokinetics of Propofol in Adult Volunteers . 
The Journal of the American Society of Anesthesiologists 88(5):p 1170-1182, May 1, 1998. 
| DOI: 10.1097/00000542-199805000-00006

https://www.researchgate.net/publication/232211950_The_Influence_of_Method_of_Administration_and_Covariates_on_the_Pharmacokinetics_of_Propofol_in_Adult_Volunteers
'''
import leanbodymass as LBM

def schnider_params_calc(a, w, h):
    v1 = 4.27 #L 
    v2 = 18.9 - 0.391 * (a - 53) #L vd_rapid_peripheral
    v3 = 238 #L vd_slow_peripheral
    clearance_met = 1.89 + ((w - 77) * 0.0456) + ((LBM.lbm_calc - 59) * -0.0681) + ((h - 177) * 0.0264)
    clearance_rapid_periph = 1.29 - 0.024 * (a - 53)
    clearance_slow_periph = 0.836 #all clearances L min-1
    
    #to make integration easier, I calculate elim constants
    k10 = clearance_met / v1
    k12 = clearance_rapid_periph / v1
    k21 = clearance_rapid_periph /v2
    k13 = clearance_slow_periph / v1
    k31 = clearance_slow_periph / v3
    
    schneider_params = [v1,v2,v3,k10,k12,k21,k13,k31]
    
    return schneider_params