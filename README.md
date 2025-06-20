# Work Counting with Process Tensors

This project calculates numerically exact quantum work statistics using process tensors.

---

## Project Overview

This code was used to calculate the results presented in (link to paper).

---

## Initial Setup

To run the code, follow these instructions:

**(1) Set up a virtual environment (recommended):**

**(2) Install requirements:**
Make sure you have Python 3.12 installed. Then, from the root of the repository, run:

```bash
pip install -r requirements.txt* 
```

**(3) Create Data Folders:** 
Navigate to the `data` folder and run the following command:
```bash
    bash setup_data_folder.sh
```

---

## Tests

To ensure everything is set up correctly, run the following tests:

* **Dynamics Test:** From within the `tests` folder, execute:
    ```bash
    bash run-dynamics-tests.sh
    ```
    This script will first create an influence functional, then run the system dynamics, and finally plot these dynamics.
* **Work Characteristic Function (WCF) Test:** From within the `tests` folder, run:
    ```bash
    bash run-wcf-tests.sh
    ```
    This will calculate the WCF using the same parameters as the influence functional and system dynamics from the previous test.
* **Verify Test Results:** To confirm your tests ran properly, compare the generated figures to the `expected-dynamics.png` and `expected-wcf.png` files located in the `tests` folder.
* **Moments Test:** TODO: Add instructions for running the moments test.

---

## List of Parameters

This project utilizes two primary categories of parameters: **physical parameters** and **numerical parameters**.

* **Physical parameters** describe the system you are simulating, such as energy scales and time.
* **Numerical parameters** control how the simulation is performed, including the time step and SVD threshold.

The goal is to achieve **convergence** of the results with respect to the numerical parameters for a given set of physical parameters. When convergence is achieved, it indicates that the results are robust and not artifacts of the simulation method.

Below is a breakdown of the physical and numerical parameters used in the bash scripts. Descriptions of other parameters used elsewhere in the project can be found in their respective files.

### Physical Parameters

These are broken down into system, bath, and shared parameters.

#### System Parameters:

* **`EPSMAX`** (float): Maximum energy splitting of the TLS.
* **`EPS0`** (float): Minimum energy splitting of the TLS.
* **`TP`** (float): Protocol time, representing the duration for which the TLS is driven.

#### Bath Parameters:

* **(Note: This project uses an underdamped Drude-Lorentz spectral density. You can modify this by changing the spectral density function and corresponding dependent functions within `src/LandauerErasure_functions.py`.)**
* **`BETA`** (float): Inverse temperature. This serves as the base unit for all other parameters.
* **`GAMMA`** (float): Width of the spectral density.
* **`W0`** (float): Peak location of the spectral density.
* **`WC`** (float): Cutoff for the spectral density.
* **`ALPHA`** (float): Coupling strength between the system and bath.

### Numerical Parameters (Simulation Parameters)

* **`STEP_SIZE`**: Trotter step size, used for the discretization of the generalized time axis.
* **`PREC`**: Threshold for Singular Value Decomposition (SVD), expressed in units of $10^{-\text{PREC}}$.

### Other Parameters

These parameters are used elsewhere in the project's code:

* **`S`** (float): Equilibration time. This is used to perform the first measurement within the TPMP and must be chosen large enough to achieve converged results.
* **`M`** (int): Counting integer, used to calculate the WCF along the chi axis in units of $M \times$ STEP_SIZE.
* **`MAX_DIFF_ORDER`** (int): Used when calculating moments through the finite difference approximation to achieve higher-order accuracy. This should be an integer and a multiple of 2.

---

## Example Usage

I've created video tutorials on my YouTube channel that explain both how the code works and how to use it. Please visit the following link to see the playlist of videos associated with this project:

(link to playlist)


## References
This project made use of the code written by Valentin Link which performs the calculations to calculate the influence functional. 
Here is a link to the paper that the work is based upon: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.132.200403
Here is a link to his github repo: https://github.com/val-link/iTEBD-TEMPO

