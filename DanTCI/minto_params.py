'''
Minto, Charles F. MB, ChB; Schnider, Thomas W. MD; Egan, Talmage D. MD; Youngs, Elizabeth MD; Lemmens, 
Harry J. M. MD; Gambus, Pedro L. MD; Billard, Valerie MD; Hoke, John F. PhD; Moore, Katherine H. P. PharmD; Hermann, 
David J. PharmD; Muir, Keith T. PhD; Mandema, Jaap W. PhD; Shafer, Steven L. MD. 
Influence of Age and Gender on the Pharmacokinetics and Pharmacodynamics of Remifentanil: I. Model Development. 
The Journal of the American Society of Anesthesiologists 86(1):p 10-23, January 1, 1997. 
| DOI: 10.1097/00000542-199701000-00004 


https://journals.lww.com/anesthesiology/fulltext/1997/01000/influence_of_age_and_gender_on_the.4.aspx
'''



def lbm_calc(g,w,h):
    if g=='m':
        lbm= 1.1*w-128*(w/h)**2
    if g=='f':
        lbm= 1.07*w-148*(w/h)**2
    return lbm

#for Minto/remi, keo varies with age linearly
#See Anesthesiology 1997;86: 24, table 1.
keo= lambda age: 0.8781651900711107 -0.007038728622707376*age

def minto_params_calc(g,a,w,h):
    lbm= lbm_calc(g,w,h)
    #See table 3 in article mentioned above
    vd_central= 5.1 - 0.0201* (a-40) + 0.072*(lbm-55) #L
    vd_rapid_peripheral= 9.82 - 0.0811* (a-40) + 0.108*(lbm-55) #L
    vd_slow_peripheral= 5.42 #L
    clearance_met= 2.6 - 0.0162*(a-40) + 0.0191*(lbm-55)
    clearance_rapid_periph= 2.05-0.0301*(a-40)
    clearance_slow_periph= 0.076-0.00113*(a-40) #all clearances L min-1
    
    #to make integration easier, I calculate elim constants
    k10= clearance_met/vd_central
    k12= clearance_rapid_periph/vd_central
    k21= clearance_rapid_periph/vd_rapid_peripheral
    k13= clearance_slow_periph/vd_central
    k31= clearance_slow_periph/vd_slow_peripheral
    ke0= keo(a)
    model= 'Minto'
    drug_name= 'Remifentanil'
    units= (r'$ng\ ml^{-1}$', r'$\mu g\ Kg^{-1}\ min^{-1}$')
    ec50= 3.0
    
    params=[drug_name, 
            model,
            units, ec50,
            vd_central,
            vd_rapid_peripheral,
            vd_slow_peripheral,
            k10,
            k12,k21,
            k13,k31,
            ke0
            ]
    return params