'''
Marsh B, White M, Morton N, Kenny GN. Pharmacokinetic model driven infusion of propofol in children. Br J Anaesth. 1991 Jul;67(1):41-8. doi: 10.1093/bja/67.1.41. PMID: 1859758.

https://pubmed.ncbi.nlm.nih.gov/1859758/
'''


v1 = 0.228 #L/kg
v2 = 0.463 # L/kg
v3 = 2893 # L/kg

k10 = 0.119 #min^-1 (k10 = K10 + K1e)
k12 = 0.112 #min^-1
k13 = 0.042 #min^-1
k21 = 0.055 #min^-1
k31 = 0.0033 #min^-1
ke0 = 0.26 #min^-1