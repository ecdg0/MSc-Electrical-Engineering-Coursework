import numpy as np
a = np.array([0.5, 1.0, 0.7])  
b = np.array([215, 270, 160]) 
P_min = 39.0
P_max = 150.0

def compute_generation(lambda_val):
    P = (lambda_val - b) / (2 * a)
    # Apply Kuhn–Tucker conditions
    return np.clip(P, P_min, P_max)

def gradient_dual_dispatch(P_D, lambda_init=8.0, alpha=0.01, tol=1e-3, max_iter=1000):
    """
    Update rule:
      lambda_new = lambda_old + alpha * (P_D - sum(P))
      
    where P = clip( (lambda - b) / (2*a), P_min, P_max ).
    """
    lambda_val = lambda_init
    for iteration in range(max_iter):
        P = compute_generation(lambda_val)
        total_P = np.sum(P)
        error = P_D - total_P
        if abs(error) < tol:
            return P, lambda_val, iteration
        # Update lambda
        lambda_val = lambda_val + alpha * error
    return P, lambda_val, max_iter

# Solve for two load cases: 320 MW and 200 MW
for load in [320, 200]:
    P_opt, lambda_opt, iters = gradient_dual_dispatch(load, lambda_init=8.0, alpha=0.5)
    print(f"\nFor Total Load = {load} MW:")
    print(f"  - P1 = {P_opt[0]:.2f} MW")
    print(f"  - P2 = {P_opt[1]:.2f} MW")
    print(f"  - P3 = {P_opt[2]:.2f} MW")
    print(f"  - Incremental Cost (λ) = {lambda_opt:.2f} $/MWh")
    print(f"  - Converged in {iters} iterations")
