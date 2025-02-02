import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import streamlit as st

# Constants for the 3-compartment model (example values)
k10 = 0.1  # Elimination rate constant (1/hr)
k12 = 0.2  # Rate constant from C1 to C2
k21 = 0.15  # Rate constant from C2 to C1
k23 = 0.05  # Rate constant from C2 to C3
k32 = 0.1  # Rate constant from C3 to C2
alpha = 0.1  # Effect-site equilibration parameter (1/hr)

# PD model constants
E_max = 100  # Maximal BIS effect
EC50 = 2  # EC50 for BIS (μg/mL)
n = 2  # Hill coefficient

# Initial conditions (Initial concentrations in each compartment)
C1_0 = 0.5  # Initial concentration in central compartment (μg/mL)
C2_0 = 0.0  # Initial concentration in peripheral compartment (μg/mL)
C3_0 = 0.0  # Initial concentration in deep compartment (μg/mL)

# Time vector for simulation (in hours)
# Default time step and duration, adjustable via Streamlit sliders
default_duration = 10  # Default time duration (hours)
default_time_step = 0.01  # Default time step

# 3-compartment PK model system of differential equations
def pk_model(C, t, k10, k12, k21, k23, k32):
    C1, C2, C3 = C
    dC1_dt = -k10 * C1 - k12 * C1 + k21 * C2
    dC2_dt = k12 * C1 - (k21 + k23) * C2 + k32 * C3
    dC3_dt = k23 * C2 - k32 * C3
    return [dC1_dt, dC2_dt, dC3_dt]

# PD model to calculate BIS based on Ce
def pd_model(Ce, E_max, EC50, n):
    return E_max / (1 + (EC50 / Ce) ** n)

# Streamlit UI components
st.title("Real-Time Propofol PK/PD Simulation")

# Sliders for time duration and speed
duration = st.slider("Simulation Duration (hours)", min_value=1, max_value=24, value=default_duration)
time_step = st.slider("Time Step (hours)", min_value=0.01, max_value=1.0, value=default_time_step, step=0.01)

# Time vector for simulation (in hours)
time = np.arange(0, duration, time_step)

# Solve the ODEs
initial_conditions = [C1_0, C2_0, C3_0]
solution = odeint(pk_model, initial_conditions, time, args=(k10, k12, k21, k23, k32))

# Calculate Ce (Effect-site concentration)
C1 = solution[:, 0]
Ce = C1 / (1 + alpha)

# Calculate BIS values
BIS_values = pd_model(Ce, E_max, EC50, n)

# Plotting Cp, Ce, and BIS
fig, axs = plt.subplots(1, 2, figsize=(12, 6))

# Plot plasma concentration (Cp) and effect-site concentration (Ce)
axs[0].plot(time, C1, label='Cp (Central Compartment)', color='b')
axs[0].plot(time, Ce, label='Ce (Effect-site)', color='r')
axs[0].set_xlabel('Time (hours)')
axs[0].set_ylabel('Concentration (μg/mL)')
axs[0].legend()
axs[0].set_title('Plasma Concentration (Cp) and Effect-site Concentration (Ce)')

# Plot BIS
axs[1].plot(time, BIS_values, label='BIS', color='g')
axs[1].set_xlabel('Time (hours)')
axs[1].set_ylabel('BIS')
axs[1].legend()
axs[1].set_title('BIS Response')

# Show the plot in the Streamlit app
st.pyplot(fig)

