import os
import pandas as pd
import numpy as np
import itertools

def read_generator_data(filename="data2.xlsx", sheet_name=0):
    """
    Reads generator data from 'data2.xlsx'.
    Expected columns:
      Unit, x^2 term, x term, constant, Minimum MW, Maximum MW, Start, Shutdown
    """
    if not os.path.isabs(filename):
        script_dir = os.path.dirname(os.path.realpath(__file__))
        filename = os.path.join(script_dir, filename)
    
    try:
        df = pd.read_excel(filename, sheet_name=sheet_name)
    except FileNotFoundError:
        raise FileNotFoundError(f"File {filename} not found.")
    except Exception as e:
        raise Exception(f"Error reading {filename}: {e}")
    
    required_cols = ["Unit", "x^2 term", "x term", "constant", "Minimum MW", "Maximum MW", "Start", "Shutdown"]
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")
    return df

def calculate_power_for_lambda(lmbda, a, b, p_min, p_max):
    """
    Computes the power output for a generator given λ, cost coefficients a and b,
to the interval [p_min, p_max].
    Formula: P = (λ - b) / (2a)
    """
    if abs(a) < 1e-14:
        p_calc = p_max if lmbda > b else p_min
    else:
        p_calc = (lmbda - b) / (2.0 * a)
    return max(p_min, min(p_calc, p_max))

def economic_dispatch_for_state(state, P_load, data_df, tolerance=0.1, max_iterations=20, lmbda_init=8.0):
    """
    Solves the economic dispatch problem for a given commitment state (list/tuple of 0/1)
    and stage load P_load.
    Only generators with u==1 are dispatched.
    Returns:
      - optimal lambda (λ)
      - list of generation outputs for all generators (0 for off units)
      - total ED cost (fuel cost only)
      
    If the sum of maximum or minimum capacities of committed units cannot meet P_load,
    returns an infinite cost.
    """
    a_list = data_df["x^2 term"].tolist()
    b_list = data_df["x term"].tolist()
    c_list = data_df["constant"].tolist()
    p_min_list = data_df["Minimum MW"].tolist()
    p_max_list = data_df["Maximum MW"].tolist()
    
    committed_indices = [i for i, u in enumerate(state) if u == 1]
    
    # Feasibility check:
    max_possible = sum(p_max_list[i] for i in committed_indices)
    min_possible = sum(p_min_list[i] for i in committed_indices)
    if max_possible < P_load or P_load < min_possible:
        return lmbda_init, [0.0]*len(state), np.inf
    
    lambdas = []
    errors = []
    lmbda = lmbda_init
    iteration = 0

    while iteration < max_iterations:
        P_vals = [0.0]*len(state)
        total_gen = 0.0
        for i in committed_indices:
            Pi = calculate_power_for_lambda(lmbda, a_list[i], b_list[i], p_min_list[i], p_max_list[i])
            P_vals[i] = Pi
            total_gen += Pi
        
        error = P_load - total_gen
        
        if abs(error) <= tolerance:
            break
        
        lambdas.append(lmbda)
        errors.append(error)
        
        if iteration < 2:
            if error > 0:
                lmbda *= 1.1
            else:
                lmbda *= 0.9
        else:
            lmbda1, e1 = lambdas[-2], errors[-2]
            lmbda2, e2 = lambdas[-1], errors[-1]
            if abs(e2 - e1) < 1e-6:
                break
            lmbda = lmbda2 - e2 * (lmbda2 - lmbda1) / (e2 - e1)
        iteration += 1

    total_ED_cost = 0.0
    for i in committed_indices:
        Pi = P_vals[i]
        cost_i = c_list[i] + b_list[i]*Pi + a_list[i]*(Pi**2)
        total_ED_cost += cost_i

    return lmbda, P_vals, total_ED_cost

def compute_transition_cost(prev_state, curr_state, data_df):
    """
    Computes the transition cost from prev_state to curr_state.
    For each generator i:
       0->1: add the Start cost.
     1->0: add the Shutdown cost.
    """
    total_trans_cost = 0.0
    for i, (prev, curr) in enumerate(zip(prev_state, curr_state)):
        if prev == 0 and curr == 1:
            total_trans_cost += data_df.iloc[i]["Start"]
        elif prev == 1 and curr == 0:
            total_trans_cost += data_df.iloc[i]["Shutdown"]
    return total_trans_cost

def state_to_str(state):
    """Converts a state tuple/list to a string for printing."""
    return ''.join(str(int(x)) for x in state)

def dynamic_programming_unit_commitment(data_filename="data2.xlsx",
                                        stage_loads=[1000, 1400, 1800, 800],
                                        tolerance=0.1,
                                        max_iterations=20,
                                        lmbda_init=8.0):
    """
    Implements dynamic programming for the multi-stage unit commitment problem.
    Returns:
      - Optimal commitment sequence (list of states for stages 0..T)
      - A list of dictionaries with detailed results for stages 1..T,
      - The sum of each stage total cost.
    """
    data_df = read_generator_data(filename=data_filename)
    num_units = len(data_df)
    stages = len(stage_loads)
    
    # Full state space (tuples of 0/1 for each unit)
    all_states = list(itertools.product([0, 1], repeat=num_units))
    
    # For stage 1, only allow states with Unit1 and Unit2 on.
    allowed_states_stage1 = [state for state in all_states if state[0] == 1 and state[1] == 1]
    
    # Filter candidate states based on feasibility for given load.
    def feasible_states(load, candidate_states):
        feasible = []
        for state in candidate_states:
            sum_max = sum(data_df.iloc[i]["Maximum MW"] for i, u in enumerate(state) if u == 1)
            sum_min = sum(data_df.iloc[i]["Minimum MW"] for i, u in enumerate(state) if u == 1)
            if sum_max >= load and load >= sum_min:
                feasible.append(state)
        return feasible

    # Stage 0: initial state (assume all off)
    initial_state = (0, 0, 0, 0)
    # DP structure: For each stage, store for each state a tuple:
    # (cumulative cost, previous state, ED_result, transition cost incurred at this stage, stage cost, ED lambda)
    dp = [{} for _ in range(stages+1)]
    dp[0][initial_state] = (0.0, None, None, 0.0, 0.0, None)
    
    for t in range(1, stages+1):
        load = stage_loads[t-1]
        dp_t = {}
        if t == 1:
            candidate_states = allowed_states_stage1
        else:
            candidate_states = all_states
        candidate_states = feasible_states(load, candidate_states)
        
        for curr_state in candidate_states:
            best_cost = np.inf
            best_prev = None
            best_ED = None
            best_trans = None
            best_stage_cost = None
            best_lambda = None
            for prev_state, (prev_cost, _, _, _, _, _) in dp[t-1].items():
                trans_cost = compute_transition_cost(prev_state, curr_state, data_df)
                lmbda, P_vals, ED_cost = economic_dispatch_for_state(curr_state, load, data_df,
                                                                       tolerance=tolerance,
                                                                       max_iterations=max_iterations,
                                                                       lmbda_init=lmbda_init)
                if np.isinf(ED_cost):
                    continue
                stage_cost = trans_cost + ED_cost
                total_cost = prev_cost + stage_cost
                if total_cost < best_cost:
                    best_cost = total_cost
                    best_prev = prev_state
                    best_ED = (P_vals, ED_cost)
                    best_trans = trans_cost
                    best_stage_cost = stage_cost
                    best_lambda = lmbda
            if best_prev is not None:
                dp_t[curr_state] = (best_cost, best_prev, best_ED, best_trans, best_stage_cost, best_lambda)
        dp[t] = dp_t
    
    if not dp[stages]:
        raise Exception("No feasible commitment found for the given stage loads.")
    
    # Select final state with minimum cumulative cost.
    final_stage = dp[stages]
    best_final_state = min(final_stage, key=lambda s: final_stage[s][0])
    best_total_cost = final_stage[best_final_state][0]
    
    # Backtrack to retrieve the optimal commitment sequence.
    optimal_states = [None]*(stages+1)
    optimal_dispatch = [None]*(stages+1)
    optimal_lambda = [None]*(stages+1)
    stage_costs = [None]*(stages+1)
    
    optimal_states[stages] = best_final_state
    optimal_dispatch[stages] = final_stage[best_final_state][2]
    optimal_lambda[stages] = final_stage[best_final_state][5]
    stage_costs[stages] = final_stage[best_final_state][4]
    
    for t in range(stages, 0, -1):
        prev_state = dp[t][optimal_states[t]][1]
        optimal_states[t-1] = prev_state
        if t-1 > 0:
            optimal_dispatch[t-1] = dp[t-1][prev_state][2]
            optimal_lambda[t-1] = dp[t-1][prev_state][5]
            stage_costs[t-1] = dp[t-1][prev_state][4]
    
    stage_results = []
    for t in range(1, stages+1):
        state = optimal_states[t]
        result = {
            "Stage": t,
            "Commitment": state,
            "ED_lambda": dp[t][state][5],
            "Dispatch": dp[t][state][2][0],
            "ED_cost": dp[t][state][2][1],
            "Transition_cost": dp[t][state][3],
            "Stage_cost": dp[t][state][4],
            "Load": stage_loads[t-1]
        }
    
        stage_results.append(result)
    
    return optimal_states, stage_results, best_total_cost

if __name__ == "__main__":
    stage_loads = [1000.0, 1400.0, 1800.0, 800.0]
    try:
        opt_states, stage_results, full_total_cost = dynamic_programming_unit_commitment(
            data_filename="data2.xlsx",
            stage_loads=stage_loads,
            tolerance=0.1,
            max_iterations=20,
            lmbda_init=8.0
        )
    except Exception as e:
        print(f"Error during DP UC computation: {e}")
        exit(1)
    
    print("Optimal Stage-by-Stage Results:")
    for res in stage_results:
        print(f"Stage {res['Stage']}:")
        print(f"  Commitment: {state_to_str(res['Commitment'])}")
        print(f"  Load: {res['Load']} MW")
        print(f"  ED λ: {res['ED_lambda']:.4f} $/MWh")
        print("  Dispatch:")
        for i, P in enumerate(res['Dispatch']):
            print(f"    Unit {i+1}: {P:.4f} MW")
        print(f"  Transition Cost: ${res['Transition_cost']:.2f}")
        print(f"  ED Cost: ${res['ED_cost']:.2f}")
        print(f"  Stage Cost: ${res['Stage_cost']:.2f}\n")
    
    print(f"Full Total Cost over 4 stages: ${full_total_cost:.2f}")
