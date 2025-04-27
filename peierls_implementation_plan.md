# Plan: Peierls Substitution for Finite Systems in Uptight

This plan outlines the steps to implement the Peierls substitution for incorporating magnetic field effects into the `uptight` code, focusing initially on finite systems (open boundary conditions) and allowing user-selectable gauges.

**Phase 1: Implementation for Finite Systems**

1.  **Information Gathering & Code Analysis:**
    *   **Goal:** Identify the precise locations in the code for modifications and understand how coordinates are handled.
    *   **Action:** Examine:
        *   `src/lib_uptight/input_data.f90`: Input parameter handling.
        *   `src/lib_uptight/TB_ham.f90`: Hamiltonian construction, hopping term calculation.
        *   `src/lib_uptight/structure.f` / `src/lib_uptight/neighbours.f90`: Coordinate and neighbour list access.
        *   `src/lib_uptight/constants.f90`: System of units and physical constants.

2.  **Input Parameter Definition:**
    *   **Goal:** Allow users to specify the magnetic field and gauge choice.
    *   **Action:** Modify `src/lib_uptight/input_data.f90` (and potentially related data structures):
        *   Add `LOGICAL :: use_magnetic_field = .FALSE.`
        *   Add `REAL, DIMENSION(3) :: magnetic_field_vector`
        *   Add `INTEGER :: gauge_choice` (e.g., 0=None, 1=Landau Z, 2=Symmetric Z). Define constants.
        *   Ensure parameters are read from the input file.

3.  **Gauge and Phase Calculation Module:**
    *   **Goal:** Create reusable functions for calculating the vector potential and the Peierls phase.
    *   **Action:** Create `magnetic_gauge.f90` or add to an existing module:
        *   Define physical constants (`e`, `hbar`, `Phi0 = h/e`) consistent with `uptight` units.
        *   Implement `FUNCTION get_vector_potential(r, B_vector, gauge_choice) RESULT(A_vector)` returning $\mathbf{A}(\mathbf{r})$. Include logic for Landau Z and Symmetric Z gauges.
        *   Implement `FUNCTION calculate_peierls_phase(R_i, R_j, B_vector, gauge_choice) RESULT(phase_factor)` returning $\exp\left(i \frac{e}{\hbar} \int_{\mathbf{R}_i}^{\mathbf{R}_j} \mathbf{A} \cdot d\mathbf{l}\right)$.
            *   Calls `get_vector_potential`.
            *   Uses midpoint approximation for the integral: $\int \approx \mathbf{A}\left(\frac{\mathbf{R}_i + \mathbf{R}_j}{2}\right) \cdot (\mathbf{R}_j - \mathbf{R}_i)$.
            *   Returns `(1.0, 0.0)` if `use_magnetic_field` is false.

4.  **Hamiltonian Modification:**
    *   **Goal:** Apply the Peierls phase factor to the hopping terms.
    *   **Action:** Modify `src/lib_uptight/TB_ham.f90`:
        *   Import necessary functions/variables.
        *   In the hopping term calculation loop:
            *   Get coordinates $\mathbf{R}_i$, $\mathbf{R}_j$.
            *   Call `phase = calculate_peierls_phase(...)`.
            *   Modify hopping: `H_ij = H_ij_original * phase`.
            *   Apply only if `use_magnetic_field` is true.
            *   Ensure Hermiticity: `H_ji = H_ji_original * CONJG(phase)`.

5.  **Testing (Finite Systems):**
    *   **Goal:** Verify correctness and gauge invariance.
    *   **Action:**
        *   Use finite system test cases (e.g., `Qdot`, `GaN_column`).
        *   Run with Landau gauge, store eigenvalues.
        *   Run with Symmetric gauge, compare eigenvalues (must match).
        *   (Optional) Compare with analytical results if possible.

6.  **Documentation:**
    *   **Goal:** Inform users how to use the feature.
    *   **Action:**
        *   Update `DEVELOPER_MANUAL.md` or other docs.
        *   Describe new input parameters and gauge options.
        *   State clearly: **Finite systems only** for this phase.

**Plan Visualization (Simplified Flow):**

```mermaid
graph TD
    A[Start: Peierls for Finite Systems] --> B{Analyze Code Structure};
    B --> C[Identify Modification Points];
    C --> D{Define Input Params};
    D --> E[Modify input_data.f90];
    E --> F{Implement Gauge/Phase Logic};
    F --> G[Create magnetic_gauge.f90];
    G --> H[Implement Core Functions];
    H --> I{Modify Hamiltonian Construction};
    I --> J[Update TB_ham.f90];
    J --> K{Test Finite Systems};
    K --> L[Run Gauge Comparison Tests];
    L --> M[Verify Gauge Invariance];
    M --> N{Document Feature};
    N --> O[Update Manual];
    O --> P[End Phase 1];

    style P fill:#f9f,stroke:#333,stroke-width:2px