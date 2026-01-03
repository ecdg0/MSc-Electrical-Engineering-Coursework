import os
import pandas as pd
import numpy as np

def read_generator_data(filename="data1.xlsx", sheet_name=0):
    """
    Reads generator data from 'data1.xlsx' located in the same directory as this script.
    Expects columns:
      Unit | x^2 term | x term | constant | Minimum MW | Maximum MW
    """
    # Build absolute path relative to this script.
    if not os.path.isabs(filename):
        script_dir = os.path.dirname(os.path.realpath(__file__))
        filename = os.path.join(script_dir, filename)
    
    try:
        df = pd.read_excel(filename, sheet_name=sheet_name)
    except FileNotFoundError:
        raise FileNotFoundError(f"File {filename} not found.")
    except Exception as e:
        raise Exception(f"Error reading {filename}: {e}")
    
    # Check that required columns exist.
    required_cols = ["Unit", "x^2 term", "x term", "constant", "Minimum MW", "Maximum MW"]
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")
    
    return df

def calculate_power_for_lambda(lmbda, a, b, p_min, p_max):
    """
    Computes the power output for a generator given λ, cost coefficients a and b,
    and then clamps the result between p_min and p_max.
    
    P = (λ - b) / (2a)
    """
    if abs(a) < 1e-14:
        # Nearly linear cost: choose max if λ > b, else min.
        p_calc = p_max if lmbda > b else p_min
    else:
        p_calc = (lmbda - b) / (2.0 * a)
    
    # Clamp to feasible limits.
    return max(p_min, min(p_calc, p_max))

def economic_dispatch_secant_method(filename="data1.xlsx",
                                    P_load=600.0,
                                    tolerance=0.1,
                                    max_iterations=20,
                                    lmbda_init=8.0):
    """
    Implements economic dispatch using an initial 10%-based adjustment followed by
    the secant method to refine the Lagrange multiplier (λ).
    Process:
      1. Compute power outputs: P_i = (λ - b_i) / (2a_i), clamped to [Pmin, Pmax].
      2. Compute error = P_load - sum(P_i).
      3. For the first two iterations, update λ by +-10% (if error is positive, increase λ; if negative, decrease λ).
      4. Then use the secant method to update λ:
             λ_new = λ2 - e2 * (λ2 - λ1) / (e2 - e1)
         until the error is within tolerance.
    """
    # Read generator data.
    df = read_generator_data(filename=filename)
    generator_data = df.to_dict("records")
    
    # Extract parameters for each generator.
    a_list = [row["x^2 term"] for row in generator_data]
    b_list = [row["x term"]    for row in generator_data]
    p_min_list = [row["Minimum MW"] for row in generator_data]
    p_max_list = [row["Maximum MW"] for row in generator_data]
    
    # Lists to store lambda and error for secant updates.
    lambdas = []
    errors = []
    
    lmbda = lmbda_init
    iteration = 0

    while iteration < max_iterations:
        # Calculate power outputs for each generator.
        P_values = []
        for i in range(len(a_list)):
            Pi = calculate_power_for_lambda(lmbda, a_list[i], b_list[i], p_min_list[i], p_max_list[i])
            P_values.append(Pi)
        
        total_generation = sum(P_values)
        error = P_load - total_generation
        
        print(f"Iteration {iteration+1}: λ = {lmbda:.4f}, Total Gen = {total_generation:.4f} MW, Error = {error:.4f} MW")
        
        # Check for convergence.
        if abs(error) <= tolerance:
            print("Convergence achieved.")
            break
        
        # Store lambda and error.
        lambdas.append(lmbda)
        errors.append(error)
        
        # For the first two iterations, adjust λ by ±10%.
        if iteration < 2:
            if error > 0:
                lmbda *= 1.1  # Increase λ if generation is low.
            else:
                lmbda *= 0.9  # Decrease λ if generation is high.
        else:
            # Use the last two points (λ1, e1) and (λ2, e2) for secant update.
            lmbda1, e1 = lambdas[-2], errors[-2]
            lmbda2, e2 = lambdas[-1], errors[-1]
            # Avoid division by zero in the secant formula.
            if abs(e2 - e1) < 1e-6:
                print("Secant method division by zero; small change in error. Exiting.")
                break
            lmbda_new = lmbda2 - e2 * (lmbda2 - lmbda1) / (e2 - e1)
            lmbda = lmbda_new
        
        iteration += 1
    
    else:
        print("Warning: Maximum iterations reached without full convergence.")
    
    print("\nFinal dispatch schedule:")
    for i, Pi in enumerate(P_values):
        print(f"  Unit {i+1}: P = {Pi:.4f} MW")
    print(f"Final λ = {lmbda:.4f} $/MWh, Total Generation = {total_generation:.4f} MW, Load = {P_load:.4f} MW, Error = {error:.4f} MW")
    
    return lmbda, P_values

if __name__ == "__main__":
    # Preview the generator data.
    try:
        data_df = read_generator_data()
        print("Preview of generator data:")
        print(data_df.head(), "\n")
    except Exception as e:
        print(f"Error reading generator data: {e}")
        exit(1)

    # Run ED
    optimal_lambda, schedule = economic_dispatch_secant_method(
        filename="data1.xlsx",
        P_load=600.0,        
        tolerance=0.1,      
        max_iterations=20,  
        lmbda_init=8.0    
    )
    # Compute total cost 
    total_cost = 0.0
    # read the generator data
    data_df = read_generator_data(filename="data1.xlsx")
    for i, P in enumerate(schedule):
        # For each generator, extract its cost coefficients.
        const = data_df.iloc[i]["constant"]
        b = data_df.iloc[i]["x term"]
        a = data_df.iloc[i]["x^2 term"]
        # Compute the cost for this generator.
        cost_i = const + b * P + a * (P ** 2)
        total_cost += cost_i

    print("Total Cost: ${:.2f}".format(total_cost))
