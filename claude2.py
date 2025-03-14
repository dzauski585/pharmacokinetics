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
    keo = 0.46  # 1/min
    
    # PD parameters (Hill function for BIS)
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
    t_start = 0
    t_end = total_time  # minutes
    t = np.arange(t_start, t_end + dt, dt)
    
    # Create infusion rate array and track bolus times
    infusion_rates = np.zeros_like(t)
    bolus_times = []
    
    # Process dosing regimen
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
        dy0dt = inf_rate * weight - k10 * y[0] - k12 * y[0] - k13 * y[0] + k21 * y[1] + k31 * y[2]
        dy1dt = k12 * y[0] - k21 * y[1]
        dy2dt = k13 * y[0] - k31 * y[2]
        dy3dt = keo * (c1 - y[3])
        
        return [dy0dt, dy1dt, dy2dt, dy3dt]
    
    # Interpolation function for infusion rates
    def infusion_func(t_val):
        idx = int(t_val / dt)
        if idx < 0 or idx >= len(infusion_rates):
            return 0
        return infusion_rates[idx]
    
    # Initial conditions
    y0 = [0, 0, 0, 0]
    
    # Solve ODE system
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
    
    # Extract results
    amounts = solution.y
    cp = amounts[0] / (v1 * 1000)  # Plasma concentration in mcg/mL (v1 is in L)
    ce = amounts[3]                # Effect site concentration already in mcg/mL
    
    # Calculate BIS using Hill equation
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

def main():
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
        {'type': 'infusion', 'time': 20.1, 'rate': 100, 'duration': 99.9},  # mcg/kg/min until 120 min
    ]
    
    # Simulation parameters
    total_time = 180  # minutes
    time_step = 0.1   # minutes
    
    # === USER EDITABLE SECTION - END ===
    
    # Get model parameters using patient demographics
    model_params = get_schneider_parameters(weight=weight, height=height, age=age, sex=sex)
    
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
    main()