def lbm_calc(g, w, h):
    if g =='m':
        lbm = 1.1 * w - 128 * (w / h)**2
    if g =='f':
        lbm = 1.07 * w - 148 * (w / h)**2
        
    return lbm