"""
1-D Transient Heat Conduction Through a Wall
------------------------------------------------

This script calculates the transient temperature distribution through
a 150 mm thick wall using an implicit finite-difference formulation.

Physical model
--------------
- One-dimensional heat conduction through the wall thickness
- Six computational nodes
- Constant thermal properties
- Constant outside and inside air temperatures
- Constant outside heat flux
- Constant convective heat-transfer coefficient

Boundary conditions
-------------------
Outside surface (Node 0):
    - Outside air temperature = 35 °C
    - Applied heat flux       = 700 W/m²
    - Convection coefficient  = 15 W/m²-K

Inside surface (Node 5):
    - Inside air temperature  = 25 °C
    - Convection coefficient  = 15 W/m²-K

Initial condition
-----------------
The entire wall initially has a temperature of 25 °C.

Numerical method
----------------
An implicit formulation is used. At each time step, the temperatures
at all six nodes at the new time step are solved simultaneously using
SymPy.
"""

import numpy as np
import sympy as sym
import matplotlib.pyplot as plt


# =============================================================================
# 1. INPUT PARAMETERS
# =============================================================================

# ---- Wall geometry ----------------------------------------------------------

wall_thickness = 0.15          # Wall thickness [m]
wall_area = 1.0                # Wall area considered [m²]

number_of_nodes = 6            # Number of temperature nodes
dx = wall_thickness / (number_of_nodes - 1)

# ---- Material properties ----------------------------------------------------

thermal_conductivity = 0.030   # Thermal conductivity [W/m-K]
density = 120                   # Density [kg/m³]
specific_heat = 70              # Specific heat capacity [J/kg-K]

# Thermal diffusivity:
# alpha = k / (rho * cp)
thermal_diffusivity = (
    thermal_conductivity / (density * specific_heat)
)

print("Thermal diffusivity =", thermal_diffusivity, "m²/s")


# ---- Boundary conditions ---------------------------------------------------

outside_temperature = 35       # Outside air temperature [°C]
inside_temperature = 25        # Inside air temperature [°C]

outside_heat_flux = 700        # Applied outside heat flux [W/m²]

convective_coefficient = 15    # Convective heat-transfer coefficient [W/m²-K]


# ---- Initial condition -----------------------------------------------------

initial_temperature = 30       # Initial wall temperature [°C]


# ---- Time settings ----------------------------------------------------------

time_step = 120                # Time step [s]
number_of_time_steps = 7       # Number of time steps to calculate

print("Time step =", time_step, "s")


# =============================================================================
# 2. NODE LOCATIONS
# =============================================================================

# Node locations through the wall thickness
node_positions = np.linspace(
    0.0,
    wall_thickness,
    number_of_nodes
)


# =============================================================================
# 3. DEFINE SYMBOLIC TEMPERATURE VARIABLES
# =============================================================================

# Temperatures at the NEW time step.
#
# t0 = temperature at outside surface
# t1 = temperature at node 1
# ...
# t5 = temperature at inside surface

t0, t1, t2, t3, t4, t5 = sym.symbols(
    "t0 t1 t2 t3 t4 t5"
)

new_temperatures = (t0, t1, t2, t3, t4, t5)


# =============================================================================
# 4. IMPLICIT SOLVER
# =============================================================================

def solve_implicit(previous_temperatures):
    """
    Calculate temperatures at the new time step.

    Parameters
    ----------
    previous_temperatures : list
        Temperatures at the previous time step:
        [T0, T1, T2, T3, T4, T5]

    Returns
    -------
    solution : dict
        Dictionary containing the calculated temperatures at the
        new time step.
    """

    # -------------------------------------------------------------------------
    # Outside boundary node (Node 0)
    #
    # Heat entering from outside by:
    #   1. Applied heat flux
    #   2. Convection from outside air
    #
    # Heat conducted from Node 0 to Node 1
    # -------------------------------------------------------------------------

    equation_0 = sym.Eq(
        outside_heat_flux
        + convective_coefficient * (outside_temperature - t0)
        + thermal_conductivity * (t1 - t0) / dx,

        thermal_conductivity * dx
        * (t0 - previous_temperatures[0])
        / (thermal_diffusivity * 2 * time_step)
    )


    # -------------------------------------------------------------------------
    # Interior node 1
    # -------------------------------------------------------------------------

    equation_1 = sym.Eq(
        thermal_conductivity * (t0 - t1) / dx
        + thermal_conductivity * (t2 - t1) / dx,

        thermal_conductivity * dx
        * (t1 - previous_temperatures[1])
        / (thermal_diffusivity * 2 * time_step)
    )


    # -------------------------------------------------------------------------
    # Interior node 2
    # -------------------------------------------------------------------------

    equation_2 = sym.Eq(
        thermal_conductivity * (t1 - t2) / dx
        + thermal_conductivity * (t3 - t2) / dx,

        thermal_conductivity * dx
        * (t2 - previous_temperatures[2])
        / (thermal_diffusivity * 2 * time_step)
    )


    # -------------------------------------------------------------------------
    # Interior node 3
    # -------------------------------------------------------------------------

    equation_3 = sym.Eq(
        thermal_conductivity * (t2 - t3) / dx
        + thermal_conductivity * (t4 - t3) / dx,

        thermal_conductivity * dx
        * (t3 - previous_temperatures[3])
        / (thermal_diffusivity * 2 * time_step)
    )


    # -------------------------------------------------------------------------
    # Interior node 4
    # -------------------------------------------------------------------------

    equation_4 = sym.Eq(
        thermal_conductivity * (t3 - t4) / dx
        + thermal_conductivity * (t5 - t4) / dx,

        thermal_conductivity * dx
        * (t4 - previous_temperatures[4])
        / (thermal_diffusivity * 2 * time_step)
    )


    # -------------------------------------------------------------------------
    # Inside boundary node (Node 5)
    #
    # Heat conducted from Node 4 to Node 5
    # Heat transferred from Node 5 to inside air by convection
    # -------------------------------------------------------------------------

    equation_5 = sym.Eq(
        thermal_conductivity * (t4 - t5) / dx
        + convective_coefficient * (inside_temperature - t5),

        thermal_conductivity * dx
        * (t5 - previous_temperatures[5])
        / (thermal_diffusivity * 2 * time_step)
    )


    # -------------------------------------------------------------------------
    # Solve the six simultaneous equations for the six unknown
    # temperatures at the new time step.
    # -------------------------------------------------------------------------

    solution = sym.solve(
        [
            equation_0,
            equation_1,
            equation_2,
            equation_3,
            equation_4,
            equation_5
        ],
        new_temperatures
    )

    return solution


# =============================================================================
# 5. CONVERT SYMPY SOLUTION TO A LIST
# =============================================================================

def get_temperature_list(solution):
    """
    Convert the SymPy solution dictionary into a list:

    [T0, T1, T2, T3, T4, T5]
    """

    return [
        float(solution[t0]),
        float(solution[t1]),
        float(solution[t2]),
        float(solution[t3]),
        float(solution[t4]),
        float(solution[t5])
    ]


# =============================================================================
# 6. INITIAL CONDITION
# =============================================================================

# All six wall nodes initially have the same temperature.
previous_temperatures = [
    initial_temperature
] * number_of_nodes


# =============================================================================
# 7. TIME-STEPPING CALCULATION
# =============================================================================

temperature_history = []

for step in range(1, number_of_time_steps + 1):

    # Solve for temperatures at the current time step
    solution = solve_implicit(previous_temperatures)

    # Convert the solution to a simple list
    current_temperatures = get_temperature_list(solution)

    # Actual simulation time
    current_time = step * time_step

    # Store the result
    temperature_history.append(
        {
            "time": current_time,
            "temperatures": current_temperatures
        }
    )

    # Display results
    print(
        f"\nTime = {current_time} s"
    )

    print(
        "Temperatures [°C] =",
        current_temperatures
    )

    # The current solution becomes the previous solution
    # for the next time step.
    previous_temperatures = current_temperatures


# =============================================================================
# 8. PLOT TEMPERATURE DISTRIBUTION THROUGH THE WALL
# =============================================================================

plt.figure(figsize=(7, 5))

for result in temperature_history:

    plt.plot(
        node_positions,
        result["temperatures"],
        marker="o",
        label=f"{result['time']} s"
    )

plt.xlabel("Distance through wall [m]")
plt.ylabel("Temperature [°C]")

plt.title("Transient Temperature Distribution Through Wall")

plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()
