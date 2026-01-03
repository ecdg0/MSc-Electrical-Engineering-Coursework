#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd


def read_excel_files():

    base_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Construct full file paths for Excel files
    trans_file = os.path.join(base_dir, "Transmission_Line_Data.xlsx")
    load_file  = os.path.join(base_dir, "Load_Data.xlsx")
    gen_file   = os.path.join(base_dir, "Generator_Data.xlsx")
    
    try:
        trans_df = pd.read_excel(trans_file)
        trans_df.rename(columns=lambda x: x.strip(), inplace=True)
    except FileNotFoundError:
        raise FileNotFoundError(f"Error: '{trans_file}' not found. "
                                "Please ensure the file exists in the script's directory "
                                "or provide the correct file path.")
    
    try:
        load_df = pd.read_excel(load_file)
        load_df.rename(columns=lambda x: x.strip(), inplace=True)
    except FileNotFoundError:
        raise FileNotFoundError(f"Error: '{load_file}' not found. "
                                "Please ensure the file exists in the script's directory "
                                "or provide the correct file path.")
    
    try:
        gen_df = pd.read_excel(gen_file)
        gen_df.rename(columns=lambda x: x.strip(), inplace=True)
    except FileNotFoundError:
        raise FileNotFoundError(f"Error: '{gen_file}' not found. "
                                "Please ensure the file exists in the script's directory "
                                "or provide the correct file path.")
    
    return trans_df, load_df, gen_df

def build_lines_list(trans_df):
    """
    From the Transmission_Line_Data Excel file, create a list of lines.
    Each line is a tuple (from_bus, to_bus, X, flow_limit) with buses converted
    to 0-index
    """
    lines = []
    for index, row in trans_df.iterrows():
        from_bus = int(row["From Bus"]) - 1  # convert to 0-index
        to_bus   = int(row["To Bus"]) - 1
        X        = float(row["X (pu)"])
        flow_lim = float(row["Flow Limit"])
        lines.append((from_bus, to_bus, X, flow_lim))
    return lines

def build_injection_vector(num_buses, load_df, gen_df, slack_bus=0):
    """
    Build the net injection vector (in MW) for each bus.
    For each generator (except the slack bus) add its setpoint.
    For each load, subtract its consumption.
    Finally, adjust the slack bus injection to balance the system.
    """
    injections = np.zeros(num_buses)
    # Add generator injections:
    for index, row in gen_df.iterrows():
        bus = int(row["Generator Bus"]) - 1
        Pset = float(row["Pset (MW)"])
        if bus != slack_bus:
            injections[bus] += Pset
    # Subtract loads:
    for index, row in load_df.iterrows():
        bus = int(row["Load Bus"]) - 1
        P_load = float(row["P_load (MW)"])
        injections[bus] -= P_load
    # Adjust slack bus injection so that total net injection is zero.
    injections[slack_bus] = -np.sum(np.delete(injections, slack_bus))
    return injections


# 2. DC Power Flow Functions

def build_B_matrix(num_buses, lines):
    """
    Build the full susceptance matrix B (size: num_buses x num_buses) from the 
    transmission line data. For a line between buses i and j with reactance X:
      - Off-diagonals: B[i,j] = -1/X, B[j,i] = -1/X.
      - Diagonals: B[i,i] = sum(1/X) for all lines incident on bus i.
    """
    B = np.zeros((num_buses, num_buses))
    for line in lines:
        i, j, X, _ = line
        susceptance = 1.0 / X
        B[i, j] -= susceptance
        B[j, i] -= susceptance
        B[i, i] += susceptance
        B[j, j] += susceptance
    return B

def solve_dc_power_flow(B, injections, slack_bus=0):
    """
    Solve the DC load flow equations: P = B * theta, with the slack bus angle fixed to 0.
    Assumes a system base of 100 MW: injections (in MW) are converted to per unit by dividing by 100.
    The solution is obtained in radians and then converted to degrees before returning.
    """
    # System base (MW)
    S_base = 100.0
    # Convert injections to per unit
    injections_pu = injections / S_base
    
    num_buses = len(injections_pu)
    indices = list(range(num_buses))
    indices.remove(slack_bus)
    B_red = B[np.ix_(indices, indices)]
    P_red = injections_pu[indices]
    theta_red = np.linalg.solve(B_red, P_red)
    
    theta_pu = np.zeros(num_buses)
    for idx, bus in enumerate(indices):
        theta_pu[bus] = theta_red[idx]
    
    # Convert from radians to degrees
    theta_deg = theta_pu * (180.0 / np.pi)
    return theta_deg


def calculate_line_flows(theta, lines, S_base=100.0):
    """
    Compute the power flow on each line using: P_ij = (theta_i - theta_j) / X.
    """
    
    flows = []
    for line in lines:
        i, j, X, _ = line
        flow_pu = (theta[i] - theta[j]) / X  # per unit flow
        flow_MW = flow_pu * S_base           # convert to MW
        flows.append(flow_MW)
    return flows

def compute_performance_index(flows, lines):
    """
    Compute a performance index (PI) as the sum of the normalized squared flows:
      PI = sum((|P_l| / FlowLimit_l)^2) over all lines.
    """
    PI = 0.0
    for flow, line in zip(flows, lines):
        _, _, _, flow_limit = line
        PI += (abs(flow) / flow_limit)**2
    return PI


# 3. Contingency Analysis Functions

def contingency_analysis(num_buses, B, injections, lines, slack_bus=0):
    """
    For each candidate line outage, remove the line, resolve the DC load flow,
    and compute the change in performance index (ΔPI). The performance index (PI)
    is defined as:
         PI = sum((|P_l| / FlowLimit_l)^2)
    over all lines.
    
    Since solve_dc_power_flow returns angles in degrees, we convert them to radians
    before computing line flows.
    
    Returns:
      - contingency_results: a sorted list of tuples (k, line, ΔPI) where k is the index of the outaged line.
      - PI_base: the performance index for the base case.
    """
    # Base case: get bus angles (in degrees) and convert to radians for flows.
    theta_base_deg = solve_dc_power_flow(B, injections, slack_bus)
    theta_base_rad = theta_base_deg * (np.pi / 180.0)
    flows_base = calculate_line_flows(theta_base_rad, lines)
    PI_base = compute_performance_index(flows_base, lines)
    
    contingency_results = []
    for k, line in enumerate(lines):
        # Remove the kth line from the list.
        lines_outage = [l for idx, l in enumerate(lines) if idx != k]
        B_outage = build_B_matrix(num_buses, lines_outage)
        try:
            theta_outage_deg = solve_dc_power_flow(B_outage, injections, slack_bus)
            theta_outage_rad = theta_outage_deg * (np.pi / 180.0)
            flows_outage = calculate_line_flows(theta_outage_rad, lines_outage)
            PI_outage = compute_performance_index(flows_outage, lines_outage)
            delta_PI = PI_outage - PI_base
        except np.linalg.LinAlgError:
            # If the outage leads to a singular network, mark it as infinitely severe.
            delta_PI = float('inf')
        contingency_results.append((k, line, delta_PI))
    contingency_results.sort(key=lambda x: x[2], reverse=True)
    return contingency_results, PI_base



# Sensitivity (PTDF and LODF) Functions

def get_reduced_mapping(num_buses, slack_bus=0):
    """
    Create a mapping from full bus indices to the reduced indices (excluding the slack bus).
    """
    mapping = {}
    reduced_index = 0
    for bus in range(num_buses):
        if bus == slack_bus:
            continue
        mapping[bus] = reduced_index
        reduced_index += 1
    return mapping

def compute_reduced_B_matrix(B, slack_bus=0):
    """
    Form the reduced B matrix by removing the row and column corresponding to the slack bus.
    """
    num_buses = B.shape[0]
    indices = list(range(num_buses))
    indices.remove(slack_bus)
    return B[np.ix_(indices, indices)]

def compute_PTDFs_for_transaction(B, slack_bus, lines, m, n):
    """
    Compute the Power Transfer Distribution Factors (PTDFs) for a given transaction
    from bus m (injection) to bus n (withdrawal).
    """
    num_buses = B.shape[0]
    mapping = get_reduced_mapping(num_buses, slack_bus)
    B_red = compute_reduced_B_matrix(B, slack_bus)
    d = np.zeros(B_red.shape[0])
    if m != slack_bus:
        d[mapping[m]] = 1
    if n != slack_bus:
        d[mapping[n]] = -1
    theta_sens = np.linalg.solve(B_red, d)
    
    PTDFs = []
    for line in lines:
        i, j, X, _ = line
        angle_i = theta_sens[mapping[i]] if i != slack_bus else 0.0
        angle_j = theta_sens[mapping[j]] if j != slack_bus else 0.0
        ptdf = (angle_i - angle_j) / X
        PTDFs.append(ptdf)
    return PTDFs

def get_incidence_vector(line, slack_bus, num_buses, mapping):
    """
    Construct the incidence vector for a given line.
    """
    i, j, _, _ = line
    n_red = num_buses - 1
    a = np.zeros(n_red)
    if i != slack_bus:
        a[mapping[i]] = 1
    if j != slack_bus:
        a[mapping[j]] = -1
    return a

def compute_LODF_matrix(B, slack_bus, lines):
    """
    Compute the LODF matrix for the network using the formula:
    
        LODF_{l,k} = (1/x_l * a_l^T X a_k) / (1 - (a_k^T X a_k)/x_k)
        
    for each candidate outage branch k and affected branch l.
    By convention, the self‐LODF (for the candidate outage branch itself) is forced to -1.
    
    Parameters:
      B         : full susceptance matrix
      slack_bus : index of the slack bus
      lines     : list of tuples (from_bus, to_bus, reactance, flow_limit)
      
    Returns:
      A square NumPy array of LODF values.
    """
    num_buses = B.shape[0]
    mapping = get_reduced_mapping(num_buses, slack_bus)
    B_red = compute_reduced_B_matrix(B, slack_bus)
    X = np.linalg.inv(B_red)
    num_lines = len(lines)
    LODF = np.zeros((num_lines, num_lines))
    
    # Precompute incidence vectors for each line in the reduced space.
    a_vectors = []
    for line in lines:
        a = get_incidence_vector(line, slack_bus, num_buses, mapping)
        a_vectors.append(a)
    
    # Loop over candidate outage branch k and affected branch l.
    for k in range(num_lines):
        a_k = a_vectors[k]
        x_k = lines[k][2]  # reactance of candidate outage branch k
        chi_k = a_k @ X @ a_k  # scalar: a_k^T X a_k
        denom = 1.0 - (chi_k / x_k)
        for l in range(num_lines):
            a_l = a_vectors[l]
            x_l = lines[l][2]
            LODF[l, k] = (1.0 / x_l) * (a_l @ X @ a_k) / denom
        # Force the self-LODF to -1 by convention.
        LODF[k, k] = -1.0
    return LODF
def analyze_line_outages_flows(B, slack_bus, lines, injections):
    """
    For each line in 'lines' (treated as a candidate outage),
    remove that line, solve the DC load flow for the reduced network,
    compute line flows, and check for limit violations.

    Returns a dictionary, where each key is a (from_bus, to_bus) tuple
    indicating the outaged line, and each value is a list of tuples describing
    limit violations. Each violation tuple has:
       (affected_line_label, flow_MW, limit, violation_ratio).
    If no lines are violated, the list is empty.
    If the network becomes disconnected, the value is the string "Disconnected".
    """

    num_buses = B.shape[0]
    outage_results = {}

    for k, candidate_line in enumerate(lines):
        from_bus_cand, to_bus_cand, _, _ = candidate_line

        # Build a new line list excluding the candidate outage line.
        lines_outage = [
            ln for idx, ln in enumerate(lines)
            if idx != k
        ]

        # Rebuild susceptance matrix B_outage
        B_outage = build_B_matrix(num_buses, lines_outage)

        try:
            # Solve DC load flow
            theta_deg = solve_dc_power_flow(B_outage, injections, slack_bus)
            theta_rad = theta_deg * (np.pi / 180.0)

            # Calculate line flows in MW (assuming e.g. S_base=100 if needed)
            flows_outage = calculate_line_flows(theta_rad, lines_outage, S_base=100.0)

            # Check each line for violations
            violations = []
            for flow_val, ln in zip(flows_outage, lines_outage):
                fbus, tbus, _, limit = ln
                if abs(flow_val) > limit:
                    ratio = abs(flow_val) / limit
                    label = f"{fbus+1}-{tbus+1}"
                    violations.append((label, flow_val, limit, ratio))

        except np.linalg.LinAlgError:
            # If the outage leads to a singular network (disconnected),
            # mark it accordingly
            violations = "Disconnected"

        # Store the results in the dictionary
        outage_label = (from_bus_cand, to_bus_cand)
        outage_results[outage_label] = violations

    return outage_results

def display_outage_flows_violations(outage_results):
    """
    Displays the dictionary returned by analyze_line_outages_flows in a readable table.
    Each row corresponds to one candidate outage. If limit violations occur,
    they are listed. Otherwise, "No violations" is shown. If disconnected, we mark it.
    """
    print("\nCONTINGENCY ANALYSIS (Line Outages and Violations)")
    print("Outage         Violated Lines (Flow MW, Limit MW, Ratio)")

    for outage, viols in outage_results.items():
        fbus, tbus = outage
        label = f"{fbus+1}-{tbus+1}"
        if viols == "Disconnected":
            viol_str = "Network Disconnected"
        elif not viols:  # empty list => no violations
            viol_str = "No violations"
        else:
            # Build a string that lists all violations
            items = []
            for (line_lbl, flow_val, limit, ratio) in viols:
                items.append(f"{line_lbl} ({flow_val:.2f}/{limit:.2f}, {ratio:.2f})")
            viol_str = "; ".join(items)
        print(f"{label:15} | {viol_str}")

# -----------------------------------------------------------------------------
# 5. Main Function
# -----------------------------------------------------------------------------
def main():
    # Read Excel files using the updated function
    trans_df, load_df, gen_df = read_excel_files()
    print("Excel files successfully loaded.")
    
    # Debug: Print column names to verify headers
    print("Transmission_Line_Data columns:", trans_df.columns.tolist())
    print("Load_Data columns:", load_df.columns.tolist())
    print("Generator_Data columns:", gen_df.columns.tolist())
    
    # Determine the total number of buses in the system.
    num_buses = int(max(trans_df["From Bus"].max(),
                        trans_df["To Bus"].max(),
                        load_df["Load Bus"].max(),
                        gen_df["Generator Bus"].max()))
    print(f"Number of buses in the network: {num_buses}")
    
    # Choose the slack bus (here, Bus 1 -> index 0)
    slack_bus = 0
    
    # Build the injection vector (in MW)
    injections = build_injection_vector(num_buses, load_df, gen_df, slack_bus)
    print("Injection vector (MW):", injections)
    
    # Build the list of transmission lines
    lines = build_lines_list(trans_df)
    
    # Build the full susceptance matrix B.
    B = build_B_matrix(num_buses, lines)
    
    # Solve the base case DC load flow.
    theta_deg = solve_dc_power_flow(B, injections, slack_bus)
    print("\nBus Voltage Angles (degrees):")
    for i, angle in enumerate(theta_deg):
        print(f"  Bus {i+1}: {angle:.2f}")
    
    # Compute line flows.
    # Convert degrees back to radians for flow calculations.
    theta_rad = theta_deg * (np.pi / 180.0)
    flows = calculate_line_flows(theta_rad, lines, S_base=100.0)
    print("\nLine Flows (MW):")
    for idx, flow in enumerate(flows):
        i, j, X, flow_limit = lines[idx]
        print(f"  Line from Bus {i+1} to Bus {j+1}: {flow:.4f}  (Limit: {flow_limit})")
    
    # Compute the base performance index (PI).
    PI_base = compute_performance_index(flows, lines)
    print(f"\nBase Performance Index (PI): {PI_base:.4f}")
    

 # Perform contingency (N-1) analysis.
    contingencies, _ = contingency_analysis(num_buses, B, injections, lines, slack_bus)
    print("\nContingency Ranking (by ΔPI):")
    for rank, (k, line, delta_PI) in enumerate(contingencies):
        i, j, X, flow_limit = line
        if delta_PI == float('inf'):
            print(f"  Outage of line {k} (Bus {i+1}-{j+1}): Network becomes disconnected")
        else:
            print(f"  Outage of line {k} (Bus {i+1}-{j+1}): ΔPI = {delta_PI:.4f}")
    
    # --- PTDF Matrix Calculation for All Unique Transactions ---
    # Generate the list of unique transactions (all unordered bus pairs m<n)
    transaction_list = []
    for m in range(num_buses):
        for n in range(m+1, num_buses):
            transaction_list.append((m, n))
    
    num_transactions = len(transaction_list)
    num_lines = len(lines)
    
    # Initialize the PTDF matrix: rows = monitored lines, columns = transactions.
    ptdf_matrix = np.zeros((num_lines, num_transactions))
    
    # For each unique transaction, compute the PTDF for each monitored line.
    for j, (m, n) in enumerate(transaction_list):
        ptdfs = compute_PTDFs_for_transaction(B, slack_bus, lines, m, n)
        for i, ptdf in enumerate(ptdfs):
            ptdf_matrix[i, j] = ptdf
    
    # Display the PTDF matrix in a formatted table.
    print("\nPOWER TRANSFER DISTRIBUTION FACTORS (PTDFs)")
    print("Affected     Transaction")
    print("Line           From(Injection) - To(Withdrawn)")
    
    # Construct and print the header row.
    header_str = "               "
    for (m, n) in transaction_list:
        header_str += f"{m+1:2d} to {n+1:2d}    "
    print(header_str)
    
    # Print each monitored line (row) with corresponding PTDF values.
    for i, line in enumerate(lines):
        from_bus, to_bus, X, flow_limit = line
        row_label = f"{from_bus+1}-{to_bus+1}"
        row_str = f"{row_label:15}"
        for j in range(num_transactions):
            row_str += f"{ptdf_matrix[i, j]:12.4f}"
        print(row_str)
    

    candidate_transactions = []
    for m in range(num_buses):
        for n in range(m+1, num_buses):
            candidate_transactions.append((m, n))
    
   # --- LODF Matrix Calculation and Display in main() ---
    
    # Compute the full square LODF matrix for the branches (each branch is a candidate outage).
    LODF_matrix = compute_LODF_matrix(B, slack_bus, lines)
    
    # Multiply by 100 to express the results as percentages.
    LODF_matrix_percent = LODF_matrix * 100.0
    
    # Build candidate outage list directly from the 'lines' list.
    candidate_transactions = []
    for line in lines:
        m, n, _, _ = line
        candidate_transactions.append((m, n))
    
    num_candidates = len(candidate_transactions)  # Should equal len(lines)
    num_monitored = len(lines)
    
    # Display the LODF matrix in a formatted table.
    print("\nLINE OUTAGE DISTRIBUTION FACTORS (LODFs) [in %]")
    print("Affected      Outage of a line")
    print("Line          From - To")
    
    # Build header row: each candidate outage branch label.
    header_str = "               "
    for (m, n) in candidate_transactions:
        header_str += f"{m+1:2d} to {n+1:2d}    "
    print(header_str)
    
    # Print each row (affected branch) with its LODF values.
    for i, line in enumerate(lines):
        from_bus, to_bus, _, _ = line
        row_label = f"{from_bus+1}-{to_bus+1}"
        row_str = f"{row_label:15}"
        for j in range(num_candidates):
            if np.isnan(LODF_matrix_percent[i, j]):
                row_str += f"{'N/A':>12}"
            else:
                row_str += f"{LODF_matrix_percent[i, j]:12.2f}"
        print(row_str)


    results = analyze_line_outages_flows(B, slack_bus, lines, injections)
    display_outage_flows_violations(results)

if __name__ == '__main__':
    main()
