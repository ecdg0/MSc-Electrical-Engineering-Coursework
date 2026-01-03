import numpy as np
import matplotlib.pyplot as plt

# Generation limits
Pmin = 20
Pmax = 125

# Incremental cost functions for Unit 1 and Unit 2
def lambda1(P1):
    return 0.8 * P1 + 160

def lambda2(P2):
    return 0.9 * P2 + 120

# Total cost functions (integrals of incremental cost functions)
# For Unit 1: F1(P1)=0.8*P1+160, so integrated cost = 0.4*P1^2 + 160*P1.
def TC1(P1):
    return 0.4 * P1**2 + 160 * P1

# For Unit 2: F2(P2)=0.9*P2+120, so integrated cost = 0.45*P2^2 + 120*P2.
def TC2(P2):
    return 0.45 * P2**2 + 120 * P2

def total_cost(P1, P2):
    return TC1(P1) + TC2(P2)

# Unconstrained solution (without limits) based on equal marginal cost:
def unconstrained_solution(L):
    # Derived from: P2 = (L+50)/2.125, and then P1 = L - P2.
    P2 = (L + 50) / 2.125
    P1 = L - P2
    lam = lambda2(P2)  # or lambda1(P1), since they are equal in the unconstrained region.
    return P1, P2, lam

# Optimal dispatch that respects generation limits.
def optimal_dispatch(L):
    # Define thresholds:
    # Lower threshold: when unconstrained P1 equals Pmin (≈82.22 MW).
    L_low = 82.2222
    # Upper threshold: when unconstrained P2 reaches Pmax (≈215.625 MW).
    L_high = 215.625
    
    if L < L_low:
        # Low load region: unconstrained solution would yield P1 < Pmin.
        P1 = Pmin
        P2 = L - Pmin
        lam = lambda2(P2)  # Marginal cost is determined by Unit 2.
    elif L <= L_high:
        # Intermediate region: unconstrained solution is feasible.
        P2 = (L + 50) / 2.125
        P1 = L - P2
        lam = lambda2(P2)  # or lambda1(P1), they are equal.
    else:
        # High load region: unconstrained solution would yield P2 > Pmax.
        P2 = Pmax
        P1 = L - Pmax
        lam = lambda1(P1)  # Marginal cost is determined by Unit 1.
    return P1, P2, lam

# --------------------------
# Plot 1: Incremental Cost Functions with Generation Limits and Example Intersection
# --------------------------
P1_range = np.linspace(Pmin, Pmax, 300)
P2_range = np.linspace(Pmin, Pmax, 300)

lambda1_vals = lambda1(P1_range)
lambda2_vals = lambda2(P2_range)

L_example = 150  # Example total load in MW
P1_ex, P2_ex, lam_ex = unconstrained_solution(L_example)

plt.figure()
plt.plot(P1_range, lambda1_vals, label='Unit 1: λ = 0.8·P₁ + 160', color='blue')
plt.plot(P2_range, lambda2_vals, label='Unit 2: λ = 0.9·P₂ + 120', color='green')
plt.axvline(Pmin, color='black', linestyle='--', label='Min Generation (20 MW)')
plt.axvline(Pmax, color='black', linestyle='--', label='Max Generation (125 MW)')
plt.axhline(lam_ex, color='red', linestyle='--', label=f'λ at L = {L_example} MW')
plt.plot(P1_ex, lambda1(P1_ex), 'bo', label=f'Intersection on Unit 1 (P₁ = {P1_ex:.2f} MW)')
plt.plot(P2_ex, lambda2(P2_ex), 'go', label=f'Intersection on Unit 2 (P₂ = {P2_ex:.2f} MW)')
plt.xlabel('Generation (MW)')
plt.ylabel('Incremental Cost ($/MWh)')
plt.title(f'Incremental Cost Functions with L = {L_example} MW')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# --------------------------
# Plot 2: Optimal Generation Allocation vs. Total Load with Load Boundaries
# --------------------------
L_vals = np.linspace(50, 250, 201)
P1_vals = []
P2_vals = []
for L in L_vals:
    p1, p2, _ = optimal_dispatch(L)
    P1_vals.append(p1)
    P2_vals.append(p2)
P1_vals = np.array(P1_vals)
P2_vals = np.array(P2_vals)

plt.figure()
plt.plot(L_vals, P1_vals, label='Unit 1 Output (P₁)', color='blue')
plt.plot(L_vals, P2_vals, label='Unit 2 Output (P₂)', color='green')
L_low = 82.2222
L_high = 215.625
plt.axvline(L_low, color='red', linestyle='--', label=f'Lower Load Bound ≈ {L_low:.2f} MW')
plt.axvline(L_high, color='red', linestyle='--', label=f'Upper Load Bound ≈ {L_high:.2f} MW')
plt.xlabel('Total Load (MW)')
plt.ylabel('Generation (MW)')
plt.title('Optimal Generation vs. Total Load')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# --------------------------
# Plot 3: Total Incremental Cost vs. Total Load
# --------------------------
TC_vals = []
for L in L_vals:
    p1, p2, _ = optimal_dispatch(L)
    TC_vals.append(total_cost(p1, p2))
TC_vals = np.array(TC_vals)

delta_L = L_vals[1] - L_vals[0]
incremental_cost = np.diff(TC_vals) / delta_L
L_mid = (L_vals[:-1] + L_vals[1:]) / 2

plt.figure()
plt.plot(L_mid, incremental_cost, 'r-', label='Total Incremental Cost')
plt.xlabel('Total Load (MW)')
plt.ylabel('Incremental Cost ($/MWh)')
plt.title('Total Incremental Cost vs. Total Load')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# --------------------------
# Plot 4: Total Production Cost vs. Total Load
# --------------------------
plt.figure()
plt.plot(L_vals, TC_vals, 'm-', label='Total Production Cost')
plt.xlabel('Total Load (MW)')
plt.ylabel('Total Production Cost ($/h)')
plt.title('Total Production Cost vs. Total Load')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
