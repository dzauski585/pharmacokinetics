def fatfreemass(g, a, w, h):
    
    bmi = w / ((h / 100)**2)
    if g == 'm':
        r = (0.88 + ((1 - 0.88)/(1 + (a / 13.4)**(-12.7)))) * (9270 * w / (6680 + 216 * bmi))
    elif g == 'f':
        r = (1.11 + ((1 - 1.11) / (1 + ((a / 7.1)**(-1.1))))) * (9270 * w / (8780 + 244 * bmi))
    
    return r 

#calculator after Al-Sallami