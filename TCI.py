import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
import leanbodymass as lbm
import schnider_params as sp

weight = 75 #kg
height = 170 #cm
age = 85 #yr
gender = 'f'
model = 'schneider'


if model == 'schneider':
    params = sp.schnider_params_calc(age, weight, height, gender)
    V1, V2, V3, k10, k12, k21, k13, k31 = params

bolus_dose_mg_per_kg = 2.5  # Bolus dose in mg/kg (adjust based on drug)

# Infusion rate (units of drug infusion in mg/h)
initial_inf_rate = 5  # Initial infusion rate in mg/h

# Pharmacodynamic parameters for effect (EC50 and Hill coefficient)
EC50 = 2  # EC50 value (mg/L)
Hill_coeff = 1.5  # Hill coefficient (steepness of the dose-response curve)

# Initial concentrations in each compartment
initial_concentrations = [0, 0, 0, 0]  # [C1, C2, C3, C4] in mg/L
V4 = 0
# Time array (in hours)
time = np.linspace(0, 24, 1000)  # 24 hours simulation

# Calculate bolus dose
bolus_dose = weight * bolus_dose_mg_per_kg  # Total bolus dose in mg
initial_concentrations = [bolus_dose / V1, 0, 0, 0]  # Initial concentrations in mg/L

# Create a list to log infusion rate changes
infusion_log = []

# Define the PK model
def pk_model(t, y, V1, V2, V3, V4, k10, k12, k13, k21, k31, initial_inf_rate, EC50, Hill_coeff):
    C1, C2, C3, C4 = y
    # Calculate the drug effect
    effect = C1**Hill_coeff / (EC50**Hill_coeff + C1**Hill_coeff)  # Hill equation
    infusion_rate = initial_inf_rate * (1 - effect)  # Adjust infusion rate

    # Log infusion rate changes
    if not infusion_log or infusion_log[-1][1] != infusion_rate:
        infusion_log.append((t, infusion_rate))

    # Differential equations
    dC1_dt = (infusion_rate / V1) - (k10 * C1) - (k12 * C1) - (k13 * C1) + (k21 * C2) + (k31 * C3)
    dC2_dt = (k12 * C1) - (k21 * C2)
    dC3_dt = (k13 * C1) - (k31 * C3)
    dC4_dt = k10 * C1  # Elimination compartment
    return [dC1_dt, dC2_dt, dC3_dt, dC4_dt]

# Solve the ODEs
solution = integrate.solve_ivp(
    pk_model, (0, 24), initial_concentrations, t_eval=time,
    args=(V1, V2, V3, V4, k10, k12, k13, k21, k31, initial_inf_rate, EC50, Hill_coeff)
)

# Extract solutions
C1, C2, C3, C4 = solution.y

# Print infusion rate log
print("Time (hours)    Infusion Rate (mg/h)")
for t, rate in infusion_log:
    print(f"{t:.3f}            {rate:.3f}")
    
# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(time, C1, label="Central Compartment (C1)", color='blue')
plt.plot(time, C2, label="Peripheral Compartment 2 (C2)", color='green')
plt.plot(time, C3, label="Peripheral Compartment 3 (C3)", color='red')
plt.plot(time, C4, label="Elimination Compartment (C4)", color='purple')
plt.title('4-Compartment PK Model with Bolus and PD Adjustment')
plt.xlabel('Time (hours)')
plt.ylabel('Concentration (mg/L)')
plt.legend(loc='best')
plt.grid(True)
plt.show()

