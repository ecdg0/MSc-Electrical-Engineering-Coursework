from docplex.mp.model import Model

def solve_opf_dc(disconnect_line=None):
    """
    DC-OPFpara red de 4 buses
    Generadores en 1,2,4; cargas en 2 (1.000 pu), 3 (1.1787 pu)
    líneas: (1,2),(1,3),(1,4),(2,3),(3,4)
    Slack:theta=0
    Inyección neta: Pk = Pgk - Pdk
    Balance
    Flujos de lineas
    Límites
    Objetivo: minimizar s_i·Pgi
    """
    # 1. Datos en pu
    buses     = [1,2,3,4]
    gens      = [1,2,4]
    Pd        = {1:0.0, 2:1.000, 3:1.1787, 4:0.0}
    s_cost    = {1:1307, 2:1211, 4:1254}
    Pmin      = {1:0.5,   2:0.375, 4:0.45}
    Pmax      = {1:2.0,   2:1.5,   4:1.8}
    all_lines = [(1,2),(1,3),(1,4),(2,3),(3,4)]
    b_sus     = 1/0.1  
    Pnm      = 0.5   
    ang_lim   = 3.1416 

    #Filtrar línea desconectada
    if disconnect_line:
        i, j = map(int, disconnect_line.split('-'))
        lines = [ℓ for ℓ in all_lines if ℓ != (i,j) and ℓ != (j,i)]
    else:
        lines = all_lines

    # 2Construcción del modelo
    mdl = Model(name='DC_OPF_with_disconnect')

    # Variables
    Pg    = {i: mdl.continuous_var(lb=Pmin[i], ub=Pmax[i], name=f'Pg{i}') for i in gens}
    theta = {i: mdl.continuous_var(lb=-ang_lim, ub=ang_lim, name=f'theta{i}') for i in buses}
    F     = {l: mdl.continuous_var(lb=-Pnm, ub=Pnm, name=f'F{l[0]}{l[1]}') for l in lines}

    # 3 Restricciones

    # Slack bus
    mdl.add_constraint(theta[1] == 0, 'slack_bus')

    # Balance nodal:
    for k in buses:
        gen_term = Pg[k] if k in Pg else 0.0
        inflow   = mdl.sum(F[(i,j)] for (i,j) in lines if j == k)
        outflow  = mdl.sum(F[(i,j)] for (i,j) in lines if i == k)
        mdl.add_constraint(gen_term - Pd[k] + inflow - outflow == 0, f'balance_{k}')

    # Ecuaciones de flujo
    for (i,j) in lines:
        mdl.add_constraint(F[(i,j)] == b_sus * (theta[i] - theta[j]), f'flow_eq_{i}_{j}')

    # 4. Objetivo: minimizar costo total
    mdl.minimize(mdl.sum(s_cost[i] * Pg[i] for i in gens))

    # 5. Resolver
    sol = mdl.solve(log_output=False)
    scenario = disconnect_line or 'ninguna'
    if not sol:
        print(f"No factible (línea desconectada: {scenario})")
        return

    # 6. Resultados
    print(f"\n--- DC-OPF (línea desconectada: {scenario}) ---")
    print(f"Costo total = {sol.objective_value:.2f} $/h")
    print("Pg (pu):", ", ".join(f"{i}={Pg[i].solution_value:.4f}" for i in gens))
    print("Flujos (pu):", ", ".join(f"{i}-{j}={F[(i,j)].solution_value:.4f}" for (i,j) in lines))
    print("Ángulos (rad):", ", ".join(f"{i}={theta[i].solution_value:.4f}" for i in buses))


if __name__ == "__main__":
    # Caso base
    solve_opf_dc()
    # desconectar línea 1-2
    solve_opf_dc(disconnect_line='1-3')
