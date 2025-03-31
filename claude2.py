import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def get_schneider_parameters(weight=70, height=170, age=40, sex='male'):
    """
    Calculate Schneider model parameters based on patient characteristics.
    
    Returns:
        dict: Dictionary containing all PK parameters
    """
    # Lean body mass calculation
    if sex.lower() == 'male':
        lbm = 1.1 * weight - 128 * (weight/height)**2
    else:  # female
        lbm = 1.07 * weight - 148 * (weight/height)**2
    
    # Volume calculations
    v1 = 4.27  # L - Central compartment
    v2 = 18.9 - 0.391 * (age - 53)  # L - Rapid peripheral compartment
    v3 = 238  # L - Slow peripheral compartment
    
    # Clearance calculations
    cl1 = 1.89 + 0.0456 * (weight - 77) - 0.0681 * (lbm - 59) + 0.0264 * (height - 177)  # L/min
    cl2 = 1.29 - 0.024 * (age - 53)  # L/min
    cl3 = 0.836  # L/min
    
    # Calculate rate constants
    k10 = cl1 / v1
    k12 = cl2 / v1
    k13 = cl3 / v1
    k21 = cl2 / v2
    k31 = cl3 / v3
    
    # Effect site parameters 
    # keo is the rate constant for equilibration between plasma and effect site. The value 0.46/min means 
    # it takes about 2-3 minutes to reach equilibrium (time to peak effect). This parameter comes from clinical 
    # studies where researchers measured both drug plasma concentrations and clinical effects (like BIS values), 
    # then calculated the time lag
    
    keo = 0.46  # 1/min
    
    # PD parameters (Hill function for BIS)
    EC50 = 2.2  # mcg/mL concentration producing 50% of maximum effect
    gamma = 1.8  # Hill coefficient steepness of concentration effect curve
    E0 = 100     # Baseline BIS value fully awake
    Emax = 0     # Maximum effect (lowest BIS) 0
    
    return {
        'v1': v1, 'v2': v2, 'v3': v3,
        'cl1': cl1, 'cl2': cl2, 'cl3': cl3,
        'k10': k10, 'k12': k12, 'k13': k13, 'k21': k21, 'k31': k31,
        'keo': keo,
        'EC50': EC50, 'gamma': gamma, 'E0': E0, 'Emax': Emax
    }
    
def get_marsh_parameters(weight=70, height=170, age=40, sex='male'):
    """Calculate Marsh model parameters based on patient characteristics."""
    # Define Marsh model parameters
    v1 = 0.228 * weight  # L - Central compartment
    v2 = 0.463 * weight  # L - Rapid peripheral compartment
    v3 = 2.893 * weight  # L - Slow peripheral compartment
    
    # Define clearances
    cl1 = 0.033 * weight  # L/min
    cl2 = 0.076 * weight  # L/min
    cl3 = 0.0064 * weight  # L/min
    
    # Calculate rate constants
    k10 = cl1 / v1
    k12 = cl2 / v1
    k13 = cl3 / v1
    k21 = cl2 / v2
    k31 = cl3 / v3
    
    # Effect site parameters and PD parameters can remain the same
    keo = 0.26  # 1/min - different for Marsh model
    EC50 = 2.2  # mcg/mL
    gamma = 1.8  # Hill coefficient
    E0 = 100     # Baseline BIS value
    Emax = 0     # Maximum effect (lowest BIS)
    
    return {
        'v1': v1, 'v2': v2, 'v3': v3,
        'cl1': cl1, 'cl2': cl2, 'cl3': cl3,
        'k10': k10, 'k12': k12, 'k13': k13, 'k21': k21, 'k31': k31,
        'keo': keo,
        'EC50': EC50, 'gamma': gamma, 'E0': E0, 'Emax': Emax
    }

def get_eleveld_parameters(weight=70, height=170, age=40, sex='male'):
    """Calculate Eleveld model parameters based on patient characteristics."""
    # Implementation of the Eleveld model parameters
    # This would include more complex calculations based on age, weight, height, sex
    # ...
    
    return {
        # Parameter dictionary similar to above
    }

def get_minto_remifentanil_parameters(weight=70, height=170, age=40, sex='male'):
    """
    Calculate Minto model parameters for remifentanil based on patient characteristics.
    
    Returns:
        dict: Dictionary containing all PK parameters for remifentanil
    """
    # Lean body mass calculation
    if sex.lower() == 'male':
        lbm = 1.1 * weight - 128 * (weight/height)**2
    else:  # female
        lbm = 1.07 * weight - 148 * (weight/height)**2
    
    # Volume calculations for remifentanil
    v1r = 5.1 - 0.0201 * (age - 40) + 0.072 * (lbm - 55)  # L - Central compartment
    v2r = 9.82 - 0.0811 * (age - 40) + 0.108 * (lbm - 55)  # L - Rapid peripheral compartment
    v3r = 5.42  # L - Slow peripheral compartment
    
    # Clearance calculations for remifentanil
    cl1r = 2.6 - 0.0162 * (age - 40) + 0.0191 * (lbm - 55)  # L/min
    cl2r = 2.05 - 0.0301 * (age - 40)  # L/min
    cl3r = 0.076 - 0.00113 * (age - 40)  # L/min
    
    # Calculate rate constants
    k10r = cl1r / v1r
    k12r = cl2r / v1r
    k13r = cl3r / v1r
    k21r = cl2r / v2r
    k31r = cl3r / v3r
    
    # Effect site parameters 
    keor = 0.595  # 1/min - faster equilibration for remifentanil
    
    # PD parameters for remifentanil
    EC50r = 19.3  # ng/mL
    gammar = 1.89  # Hill coefficient
    
    return {
        'v1r': v1r, 'v2r': v2r, 'v3r': v3r,
        'cl1r': cl1r, 'cl2r': cl2r, 'cl3r': cl3r,
        'k10r': k10r, 'k12r': k12r, 'k13r': k13r, 'k21r': k21r, 'k31r': k31r,
        'keor': keor,
        'EC50r': EC50r, 'gammar': gammar
    }
    
def get_minto_combined_parameters(weight=70, height=170, age=40, sex='male'):
    """
    Calculate parameters for combined propofol-remifentanil Minto model.
    
    Returns:
        dict: Combined dictionary with parameters for both drugs and interaction parameters
    """
    # Get propofol parameters using Schneider model
    propofol_params = get_schneider_parameters(weight, height, age, sex)
    
    # Get remifentanil parameters
    remifentanil_params = get_minto_remifentanil_parameters(weight, height, age, sex)
    
    # Interaction parameters for BIS response
    interaction_params = {
        'alpha': 1.33,    # Interaction coefficient (>1 indicates synergy)
        'E0': 100,        # Baseline BIS value
        'Emax': 0,        # Maximum effect (lowest BIS)
        'gamma': 1.43     # Steepness of the response surface
    }
    
    # Combine all parameters
    combined_params = {**propofol_params, **remifentanil_params, **interaction_params}
    
    return combined_params

def get_pk_model_function(model_name):
    """
    Return the appropriate PK model function based on model name.
    
    Args:
        model_name (str): Name of the PK model ('schneider', 'marsh', 'eleveld', etc.)
        
    Returns:
        function: The corresponding parameter calculation function
    """
    model_functions = {
        'schneider': get_schneider_parameters,
        'marsh': get_marsh_parameters,
        'eleveld': get_eleveld_parameters,
        # Add more models as needed
    }
    
    if model_name.lower() not in model_functions:
        raise ValueError(f"Unknown model: {model_name}. Available models: {', '.join(model_functions.keys())}")
    
    return model_functions[model_name.lower()]

def run_simulation(model_params, weight=70, dosing_regimen=None, total_time=180, dt=0.1):
    """
    Run PK/PD simulation using provided model parameters and dosing regimen.
    
    Args:
        model_params: Dictionary containing all required PK/PD parameters
        weight: Patient weight in kg
        dosing_regimen: List of dictionaries describing dosing events
        total_time: Total simulation time in minutes
        dt: Time step for simulation in minutes
        
    Returns:
        tuple: (time, infusion_rates, bolus_times, cp, ce, bis, wake_time)
    """
    # Extract parameters
    v1 = model_params['v1']
    k10 = model_params['k10']
    k12 = model_params['k12']
    k13 = model_params['k13']
    k21 = model_params['k21']
    k31 = model_params['k31']
    keo = model_params['keo']
    EC50 = model_params['EC50']
    gamma = model_params['gamma']
    E0 = model_params['E0']
    Emax = model_params['Emax']
    
    # Default dosing regimen if none provided
    if dosing_regimen is None:
        dosing_regimen = [
            {'type': 'bolus', 'time': 0, 'dose': 2.0, 'duration': 0.1},  # mg/kg bolus at t=0
            {'type': 'infusion', 'time': 0.1, 'rate': 166.67, 'duration': 10},  # mcg/kg/min for 10 min
            {'type': 'infusion', 'time': 10.1, 'rate': 133.33, 'duration': 10},  # mcg/kg/min for 10 min
            {'type': 'infusion', 'time': 20.1, 'rate': 100, 'duration': 99.9},  # mcg/kg/min until 120 min
        ]
    
    # Simulation time settings
    # Creates a time array from 0 to total simulation time with step size dt
    # This array becomes the x-axis for all plots and calculations
    t_start = 0
    t_end = total_time  # minutes
    t = np.arange(t_start, t_end + dt, dt)
    
    # Create infusion rate array and track bolus times
    # Creates an array to store infusion rates at each time point
    # Initializes an empty list to track bolus administration times for plotting
    infusion_rates = np.zeros_like(t)
    bolus_times = []
    
    # Process dosing regimen
    # Processes bolus doses by converting them to very short, high-rate infusions
    # For calculation purposes, it converts the bolus dose (mg/kg) to an equivalent infusion rate (mcg/kg/min)
    # Stores the calculated rates in the infusion_rates array at the appropriate time indices
    for dose in dosing_regimen:
        if dose['type'] == 'bolus':
            # Record bolus time for plotting
            bolus_times.append((dose['time'], dose['dose']))
            
            # For simulation purposes, model bolus as very short infusion
            start_idx = int(dose['time'] / dt)
            end_idx = int((dose['time'] + dose['duration']) / dt)
            if start_idx < len(infusion_rates) and end_idx < len(infusion_rates):
                # Convert mg/kg bolus over duration to mcg/kg/min
                bolus_rate = dose['dose'] * 1000 / dose['duration']  # mcg/kg/min
                infusion_rates[start_idx:end_idx+1] = bolus_rate
                
        elif dose['type'] == 'infusion':
            start_idx = int(dose['time'] / dt)
            end_idx = int((dose['time'] + dose['duration']) / dt)
            if start_idx < len(infusion_rates) and end_idx < len(infusion_rates):
                # Rate is already in mcg/kg/min
                infusion_rates[start_idx:end_idx+1] = dose['rate']
    
    # System of differential equations
    def pk_model(t, y, infusion_func):
        # y[0]: amount in central compartment (mcg)
        # y[1]: amount in peripheral compartment 1 (mcg)
        # y[2]: amount in peripheral compartment 2 (mcg)
        # y[3]: effect site concentration (mcg/mL)
        
        # Get infusion rate at current time
        inf_rate = infusion_func(t)
        
        # Compute compartment concentrations in mcg/mL
        c1 = y[0] / (v1 * 1000)  # central compartment concentration in mcg/mL (v1 is in L)
        
        # Differential equations
        # dy0dt: Rate of change in central compartment (input minus elimination and distribution, plus redistribution)
        # dy1dt: Rate of change in rapid peripheral compartment
        # dy2dt: Rate of change in slow peripheral compartment
        # dy3dt: Rate of change in effect site concentration
        # Key part: dy3dt = keo * (c1 - y[3]) models effect site equilibration based on the concentration gradient
        
        dy0dt = inf_rate * weight - k10 * y[0] - k12 * y[0] - k13 * y[0] + k21 * y[1] + k31 * y[2]
        dy1dt = k12 * y[0] - k21 * y[1]
        dy2dt = k13 * y[0] - k31 * y[2]
        dy3dt = keo * (c1 - y[3])
        
        return [dy0dt, dy1dt, dy2dt, dy3dt]
    
    # Interpolation function for infusion rates
    # Helper function that provides the infusion rate at any time point
    # Returns 0 if time is outside the simulation range
    
    def infusion_func(t_val):
        idx = int(t_val / dt)
        if idx < 0 or idx >= len(infusion_rates):
            return 0
        return infusion_rates[idx]
    
    # Initial conditions
    y0 = [0, 0, 0, 0]
    
    # Solve ODE system
    # The code uses solve_ivp with the 'BDF' (Backward Differentiation Formula) method:
    # solution = solve_ivp(    pk_model,     [t_start, t_end],     y0,     t_eval=t,     args=(infusion_func,),  
    #                          method='BDF',    rtol=1e-4,    atol=1e-6)
    # BDF is specifically chosen because:# This is a "stiff" system with processes occurring at very different time scales
    # Elimination rate (k10) may be much slower than distribution rates (k12, k13)
    # BDF handles stiff systems more efficiently than explicit methods
    # The tolerance parameters (rtol, atol) control solution accuracy
    solution = solve_ivp(
        pk_model, 
        [t_start, t_end], 
        y0, 
        t_eval=t, 
        args=(infusion_func,),
        method='BDF',
        rtol=1e-4,
        atol=1e-6
    )
    
    # Extract results and sotre in array same size as time array
    amounts = solution.y
    cp = amounts[0] / (v1 * 1000)  # Plasma concentration in mcg/mL (v1 is in L)
    ce = amounts[3]                # Effect site concentration already in mcg/mL
    
    # Calculate BIS using Hill equation
    # This is the Hill equation, a core component of the model
    # It calculates the BIS value (measure of sedation level) based on effect site concentration
    # The Hill equation is a sigmoid function commonly used in pharmacology to model drug effects
    # Origin: Named after A.V. Hill's work in the early 1900s on oxygen binding to hemoglobin
    # Mathematical form: E = Emax × C^γ / (EC50^γ + C^γ)
    bis = E0 - (E0 - Emax) * (ce**gamma) / (EC50**gamma + ce**gamma)
    
    # Find wake-up time (Ce < 1 mcg/mL or BIS > 80) after end of infusion
    wake_time = None
    
    # Find the end of the last infusion
    end_of_infusion = 0
    for dose in dosing_regimen:
        if dose['type'] == 'infusion':
            end_time = dose['time'] + dose['duration']
            if end_time > end_of_infusion:
                end_of_infusion = end_time
    
    # The code defines wake-up based on clinical thresholds:# 
    # Ce < 1 mcg/mL: Effect site concentration below awakening threshold# 
    # BIS > 80: BIS value indicating light sedation/approaching consciousness# 
    # It first identifies when the last infusion ends, then looks for the first time point after that where either condition is met
    # Find wake-up time after end of infusion
    for i in range(len(t)):
        if t[i] > end_of_infusion and (ce[i] < 1 or bis[i] > 80):
            wake_time = t[i]
            break
    
    # Create infusion display for graph (excluding boluses)
    infusion_display = np.zeros_like(infusion_rates)
    for dose in dosing_regimen:
        if dose['type'] == 'infusion':
            start_idx = int(dose['time'] / dt)
            end_idx = int((dose['time'] + dose['duration']) / dt)
            if start_idx < len(infusion_display) and end_idx < len(infusion_display):
                infusion_display[start_idx:end_idx+1] = dose['rate']
    
    return t, infusion_display, bolus_times, cp, ce, bis, wake_time
'''Combined simulation function with many similarities
def run_combined_simulation(model_params, weight=70, propofol_dosing=None, 
                           remifentanil_rate=0.05, total_time=180, dt=0.1):
    """
    Run PK/PD simulation for combined propofol-remifentanil using provided parameters.
    
    Args:
        model_params: Dictionary containing all required PK/PD parameters for both drugs
        weight: Patient weight in kg
        propofol_dosing: List of dictionaries describing propofol dosing events
        remifentanil_rate: Constant remifentanil infusion rate in mcg/kg/min
        total_time: Total simulation time in minutes
        dt: Time step for simulation in minutes
        
    Returns:
        tuple: (time, propofol_infusion_rates, bolus_times, cp_propofol, ce_propofol, 
                cp_remi, ce_remi, bis, wake_time)
    """
    # Extract propofol parameters
    v1 = model_params['v1']
    k10 = model_params['k10']
    k12 = model_params['k12']
    k13 = model_params['k13']
    k21 = model_params['k21']
    k31 = model_params['k31']
    keo = model_params['keo']
    
    # Extract remifentanil parameters
    v1r = model_params['v1r']
    k10r = model_params['k10r']
    k12r = model_params['k12r']
    k13r = model_params['k13r']
    k21r = model_params['k21r']
    k31r = model_params['k31r']
    keor = model_params['keor']
    
    # Extract interaction parameters for BIS response
    alpha = model_params['alpha']
    E0 = model_params['E0']
    Emax = model_params['Emax']
    gamma = model_params['gamma']
    EC50p = model_params['EC50']    # EC50 for propofol
    EC50r = model_params['EC50r']   # EC50 for remifentanil
    
    # Default propofol dosing regimen if none provided
    if propofol_dosing is None:
        propofol_dosing = [
            {'type': 'bolus', 'time': 0, 'dose': 2.0, 'duration': 0.1},  # mg/kg bolus at t=0
            {'type': 'infusion', 'time': 0.1, 'rate': 166.67, 'duration': 10},  # mcg/kg/min for 10 min
            {'type': 'infusion', 'time': 10.1, 'rate': 133.33, 'duration': 10},  # mcg/kg/min for 10 min
            {'type': 'infusion', 'time': 20.1, 'rate': 100, 'duration': 99.9},  # mcg/kg/min until 120 min
        ]
    
    # Simulation time settings
    t_start = 0
    t_end = total_time  # minutes
    t = np.arange(t_start, t_end + dt, dt)
    
    # Create propofol infusion rate array and track bolus times
    propofol_rates = np.zeros_like(t)
    bolus_times = []
    
    # Process propofol dosing regimen (same as original code)
    for dose in propofol_dosing:
        if dose['type'] == 'bolus':
            bolus_times.append((dose['time'], dose['dose']))
            start_idx = int(dose['time'] / dt)
            end_idx = int((dose['time'] + dose['duration']) / dt)
            if start_idx < len(propofol_rates) and end_idx < len(propofol_rates):
                bolus_rate = dose['dose'] * 1000 / dose['duration']  # mcg/kg/min
                propofol_rates[start_idx:end_idx+1] = bolus_rate
        elif dose['type'] == 'infusion':
            start_idx = int(dose['time'] / dt)
            end_idx = int((dose['time'] + dose['duration']) / dt)
            if start_idx < len(propofol_rates) and end_idx < len(propofol_rates):
                propofol_rates[start_idx:end_idx+1] = dose['rate']
    
    # Create constant remifentanil infusion rate
    remifentanil_rates = np.ones_like(t) * remifentanil_rate
    
    # System of differential equations for both drugs
    def pkpd_model(t, y, propofol_func, remifentanil_func):
        # y[0:3]: propofol amounts in compartments (mcg)
        # y[4]: propofol effect site concentration (mcg/mL)
        # y[5:7]: remifentanil amounts in compartments (mcg)
        # y[8]: remifentanil effect site concentration (ng/mL)
        
        # Get infusion rates at current time
        propofol_rate = propofol_func(t)
        remifentanil_rate = remifentanil_func(t)
        
        # Compute compartment concentrations
        cp = y[0] / (v1 * 1000)  # propofol central concentration in mcg/mL
        cr = y[4] / (v1r * 1000) * 1000  # remifentanil central concentration in ng/mL
        
        # Propofol differential equations
        dy0dt = propofol_rate * weight - k10 * y[0] - k12 * y[0] - k13 * y[0] + k21 * y[1] + k31 * y[2]
        dy1dt = k12 * y[0] - k21 * y[1]
        dy2dt = k13 * y[0] - k31 * y[2]
        dy3dt = keo * (cp - y[3])
        
        # Remifentanil differential equations (note: remifentanil rate in mcg/kg/min)
        dy4dt = remifentanil_rate * weight - k10r * y[4] - k12r * y[4] - k13r * y[4] + k21r * y[5] + k31r * y[6]
        dy5dt = k12r * y[4] - k21r * y[5]
        dy6dt = k13r * y[4] - k31r * y[6]
        dy7dt = keor * (cr - y[7])
        
        return [dy0dt, dy1dt, dy2dt, dy3dt, dy4dt, dy5dt, dy6dt, dy7dt]
    
    # Interpolation functions for infusion rates
    def propofol_func(t_val):
        idx = int(t_val / dt)
        if idx < 0 or idx >= len(propofol_rates):
            return 0
        return propofol_rates[idx]
    
    def remifentanil_func(t_val):
        idx = int(t_val / dt)
        if idx < 0 or idx >= len(remifentanil_rates):
            return 0
        return remifentanil_rates[idx]
    
    # Initial conditions for all compartments (propofol and remifentanil)
    y0 = [0, 0, 0, 0, 0, 0, 0, 0]
    
    # Solve ODE system
    solution = solve_ivp(
        pkpd_model, 
        [t_start, t_end], 
        y0, 
        t_eval=t, 
        args=(propofol_func, remifentanil_func),
        method='BDF',
        rtol=1e-4,
        atol=1e-6
    )
    
    # Extract results
    amounts = solution.y
    cp_propofol = amounts[0] / (v1 * 1000)  # Propofol plasma concentration in mcg/mL
    ce_propofol = amounts[3]                # Propofol effect site concentration in mcg/mL
    cp_remi = amounts[4] / (v1r * 1000) * 1000  # Remifentanil plasma concentration in ng/mL
    ce_remi = amounts[7]                    # Remifentanil effect site concentration in ng/mL
    
    # Calculate BIS using modified Hill equation for drug interaction
    # Normalized propofol and remifentanil concentrations
    u_propofol = ce_propofol / EC50p
    u_remi = ce_remi / EC50r
    
    # Combined drug effect using Minto model (response surface)
    U = u_propofol + u_remi - alpha * u_propofol * u_remi/(u_propofol + u_remi)
    bis = E0 - (E0 - Emax) * (U**gamma) / (1 + U**gamma)
    
    # Determine wake-up time as before, but now based on combined effect
    wake_time = None
    end_of_infusion = 0
    for dose in propofol_dosing:
        if dose['type'] == 'infusion':
            end_time = dose['time'] + dose['duration']
            if end_time > end_of_infusion:
                end_of_infusion = end_time
    
    for i in range(len(t)):
        if t[i] > end_of_infusion and bis[i] > 80:
            wake_time = t[i]
            break
    
    # Create infusion display for graph (excluding boluses)
    propofol_display = np.zeros_like(propofol_rates)
    for dose in propofol_dosing:
        if dose['type'] == 'infusion':
            start_idx = int(dose['time'] / dt)
            end_idx = int((dose['time'] + dose['duration']) / dt)
            if start_idx < len(propofol_display) and end_idx < len(propofol_display):
                propofol_display[start_idx:end_idx+1] = dose['rate']
    
    return (t, propofol_display, bolus_times, cp_propofol, ce_propofol, 
            cp_remi, ce_remi, bis, wake_time)
'''
def plot_results(t, infusion_rates, bolus_times, cp, ce, bis, wake_time=None):
    """
    Plot the results of the PK/PD simulation.
    """
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 14), sharex=True)
    
    # Plot infusion rates (excluding boluses)
    ax1.plot(t, infusion_rates, 'k-', label='Infusion Rate')
    
    # Mark bolus doses with vertical lines and annotations
    for bolus_time, bolus_dose in bolus_times:
        ax1.axvline(x=bolus_time, color='r', linestyle='--')
        ax1.annotate(f'Bolus: {bolus_dose} mg/kg', 
                    xy=(bolus_time, max(infusion_rates)/2),
                    xytext=(bolus_time + 5, max(infusion_rates)/2),
                    arrowprops=dict(facecolor='red', shrink=0.05),
                    fontsize=10)
    
    ax1.set_ylabel('Infusion Rate\n(mcg/kg/min)')
    ax1.set_title('Propofol PK/PD Simulation')
    ax1.legend()
    ax1.grid(True)
    
    # Plot concentrations
    ax2.plot(t, cp, 'b-', label='Plasma Conc. (Cp)')
    ax2.plot(t, ce, 'r-', label='Effect Site Conc. (Ce)')
    ax2.axhline(y=1, color='g', linestyle='--', label='Ce = 1 mcg/mL')
    
    if wake_time:
        ax2.axvline(x=wake_time, color='m', linestyle='--')
        ax2.annotate(f'Wake-up: {wake_time:.1f} min', 
                    xy=(wake_time, max(ce)/2),
                    xytext=(wake_time + 5, max(ce)/2),
                    arrowprops=dict(facecolor='magenta', shrink=0.05),
                    fontsize=10)
    
    ax2.set_ylabel('Concentration\n(mcg/mL)')
    ax2.set_ylim(0,8 )#ignores bolus Ce
    ax2.legend()
    ax2.grid(True)
    
    # Plot BIS
    ax3.plot(t, bis, 'g-', label='BIS')
    ax3.axhline(y=80, color='r', linestyle='--', label='BIS = 80')
    ax3.set_ylabel('BIS')
    ax3.set_xlabel('Time (min)')
    ax3.legend()
    ax3.grid(True)
    ax3.set_xlim(0, max(t))
    
    plt.tight_layout()
    return fig
''' model comparison and graphing functions
def compare_models(weight, height, age, sex, dosing_regimen=None):
    """Run simulations with different PK models and compare results."""
    models = ['schneider', 'marsh', 'eleveld']
    results = {}
    
    for model in models:
        model_function = get_pk_model_function(model)
        model_params = model_function(weight=weight, height=height, age=age, sex=sex)
        results[model] = run_simulation(model_params, weight, dosing_regimen)
    
    # Plot comparative results
    plot_model_comparison(results)
    
def plot_model_comparison(results):
    """Plot comparative results from different PK models."""
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 14), sharex=True)
    
    # Different colors for different models
    colors = {'schneider': 'b', 'marsh': 'r', 'eleveld': 'g'}
    
    # Plot Cp for each model
    for model, (t, _, _, cp, ce, bis, _) in results.items():
        ax1.plot(t, cp, color=colors[model], label=f'{model.capitalize()} Cp')
    
    # Similar plotting for Ce and BIS
    # ...
    
    plt.tight_layout()
    plt.show()
'''
'''combined model ploting function

def plot_combined_results(t, propofol_rates, bolus_times, cp_propofol, ce_propofol, 
                         cp_remi, ce_remi, bis, wake_time=None):
    """
    Plot the results of the combined propofol-remifentanil PK/PD simulation.
    """
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(12, 16), sharex=True)
    
    # Plot propofol infusion rates (same as original)
    ax1.plot(t, propofol_rates, 'k-', label='Propofol Rate')
    for bolus_time, bolus_dose in bolus_times:
        ax1.axvline(x=bolus_time, color='r', linestyle='--')
        ax1.annotate(f'Bolus: {bolus_dose} mg/kg', 
                    xy=(bolus_time, propofol_rates.max() * 0.9),
                    xytext=(bolus_time + 2, propofol_rates.max() * 0.9),
                    arrowprops=dict(arrowstyle='->'))
    ax1.set_ylabel('Propofol Rate (mcg/kg/min)')
    ax1.set_title('Propofol Infusion Rate')
    ax1.grid(True)
    
    # Plot remifentanil concentration
    ax2.plot(t, cp_remi, 'b-', label='Plasma')
    ax2.plot(t, ce_remi, 'r-', label='Effect Site')
    ax2.set_ylabel('Remifentanil Conc. (ng/mL)')
    ax2.set_title('Remifentanil Concentration')
    ax2.legend()
    ax2.grid(True)
    
    # Plot propofol concentration
    ax3.plot(t, cp_propofol, 'b-', label='Plasma')
    ax3.plot(t, ce_propofol, 'r-', label='Effect Site')
    ax3.set_ylabel('Propofol Conc. (mcg/mL)')
    ax3.set_title('Propofol Concentration')
    ax3.legend()
    ax3.grid(True)
    
    # Plot BIS
    ax4.plot(t, bis, 'g-')
    ax4.axhline(y=60, color='k', linestyle='--', label='BIS 60')
    ax4.axhline(y=40, color='k', linestyle=':', label='BIS 40')
    
    # Mark wake-up time if available
    if wake_time is not None:
        ax4.axvline(x=wake_time, color='m', linestyle='-.')
        ax4.annotate(f'Wake-up: {wake_time:.1f} min', 
                    xy=(wake_time, 60),
                    xytext=(wake_time + 2, 60),
                    arrowprops=dict(arrowstyle='->'))
    
    ax4.set_ylabel('BIS Value')
    ax4.set_xlabel('Time (minutes)')
    ax4.set_title('BIS Response')
    ax4.grid(True)
    ax4.set_ylim(0, 100)
    
    plt.tight_layout()
    plt.show()
    '''
def main(model_name='schneider'):
    # === USER EDITABLE SECTION - START ===
    
    # Patient parameters
    weight = 70      # kg
    height = 170     # cm
    age = 40         # years
    sex = 'male'     # 'male' or 'female'
    
    # Bolus and infusion parameters (already in mcg/kg/min for infusions)
    dosing_regimen = [
        {'type': 'bolus', 'time': 0, 'dose': 2.0, 'duration': 0.1},  # mg/kg bolus at t=0
        {'type': 'infusion', 'time': 0.1, 'rate': 166.67, 'duration': 10},  # mcg/kg/min for 10 min
        {'type': 'infusion', 'time': 10.1, 'rate': 133.33, 'duration': 10},  # mcg/kg/min for 10 min
        {'type': 'infusion', 'time': 20.1, 'rate': 100, 'duration': 119.9},  # mcg/kg/min until 120 min
    ]
    
    # Simulation parameters
    total_time = 180  # minutes
    time_step = 0.1   # minutes
    
    # === USER EDITABLE SECTION - END ===
    
    # Get the appropriate model function
    model_function = get_pk_model_function(model_name)
    
    # Get model parameters using patient demographics
    model_params = model_function(weight=weight, height=height, age=age, sex=sex)
    
    # Run simulation
    t, infusion_display, bolus_times, cp, ce, bis, wake_time = run_simulation(
        model_params,
        weight=weight,  # Only specify weight once
        dosing_regimen=dosing_regimen,
        total_time=total_time,
        dt=time_step
    )
    
    # Plot results
    fig = plot_results(t, infusion_display, bolus_times, cp, ce, bis, wake_time)
    plt.show()
    
    # Print wake-up time and summary
    print("\n===== SIMULATION SUMMARY =====")
    print(f"PK/PD Model: {model_name.capitalize()}")
    print(f"Patient: {weight}kg, {height}cm, {age}y, {sex}")
    
    if wake_time:
        print(f"Wake-up time: {wake_time:.2f} minutes after start")
        print(f"Recovery time: {wake_time - 120:.2f} minutes after end of infusion")
    else:
        print("No wake-up within simulation time")
    
    # Print max concentrations
    print(f"Maximum Cp: {max(cp):.2f} mcg/mL")
    print(f"Maximum Ce: {max(ce):.2f} mcg/mL")
    print(f"Minimum BIS: {min(bis):.2f}")
    
    # Print key metrics for clinicians
    print("\n--- Key Clinical Metrics ---")
    
    # Time to BIS < 60 (surgical anesthesia)
    time_to_surgical = None
    for i in range(len(t)):
        if bis[i] < 60:
            time_to_surgical = t[i]
            break
    
    # Time to Ce > 3 mcg/mL
    time_to_therapeutic = None
    for i in range(len(t)):
        if ce[i] > 3:
            time_to_therapeutic = t[i]
            break
    
    if time_to_surgical:
        print(f"Time to surgical anesthesia (BIS < 60): {time_to_surgical:.2f} minutes")
    else:
        print("Surgical anesthesia not achieved in simulation time")
    
    if time_to_therapeutic:
        print(f"Time to therapeutic Ce (> 3 mcg/mL): {time_to_therapeutic:.2f} minutes")
    else:
        print("Therapeutic Ce not achieved in simulation time")
    
    # Print drug usage
    total_drug = 0
    for dose in dosing_regimen:
        if dose['type'] == 'bolus':
            # Convert mg/kg to mg total
            total_drug += dose['dose'] * weight
        elif dose['type'] == 'infusion':
            # Convert mcg/kg/min to mg total
            # dose['rate'] * weight * duration / 1000
            total_drug += dose['rate'] * weight * dose['duration'] / 1000
    
    print(f"Total propofol used: {total_drug:.2f} mg")
    print(f"Propofol usage per minute: {total_drug/120:.2f} mg/min during maintenance")


if __name__ == "__main__":
    import argparse
    ## usage is --model 'model'
    
    parser = argparse.ArgumentParser(description='Propofol PK/PD Simulation')
    parser.add_argument('--model', type=str, default='schneider', 
                        choices=['schneider', 'marsh', 'eleveld'],
                        help='PK/PD model to use for simulation')
    # Add more arguments as needed
    
    args = parser.parse_args()
    main(model_name=args.model)