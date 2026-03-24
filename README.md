# Atmospheric Re-entry Guidance using MPC and EKF

This project implements a **nonlinear predictive guidance framework** for atmospheric re-entry using **Model Predictive Control (MPC)** combined with **Extended Kalman Filtering (EKF)** for state estimation under uncertainty.

The system models a lifting re-entry vehicle in a rotating Earth environment and evaluates performance under stochastic disturbances.

---

## 📌 Features

* 3-DOF nonlinear re-entry dynamics (ECEF frame)
* MPC-based bank angle guidance
* EKF state estimation with augmented states:

  * Ballistic coefficient scaling
  * Atmospheric density bias
  * Wind components
* Constraint handling (heat rate, dynamic pressure, g-load, bank limits)
* Reachability and footprint analysis
* Monte Carlo robustness evaluation (100 runs)

---

## 🧠 Methodology

* Nonlinear optimal control formulated using MPC
* Guidance and estimation are tightly coupled
* State uncertainty propagated via EKF covariance
* Robustness evaluated under:

  * ±20% ballistic coefficient variation
  * ±15% density bias
  * Wind disturbances (±50 m/s)
  * Navigation noise

---

## 📊 Results

The simulation framework generates:

* Ground track trajectories
* Altitude–velocity profiles
* Heat rate and g-load histories
* Bank angle profiles
* Estimation covariance evolution
* Landing footprint (3σ dispersion)
* Constraint violation statistics

---

## 📁 Project Structure

* `dynamics/`       : Vehicle dynamics models
* `guidance/`       : MPC controller
* `estimation/`     : EKF implementation
* `simulation/`     : Monte Carlo simulation routines
* `utils/`          : Helper functions
* `Figures/`        : Generated plots

---

## 🚀 How to Run

### Run Monte Carlo simulation

```matlab
run_simulation

The number of runs for MC is given in the variable 'Nsim' 
```

---

## ⚙️ Requirements

* MATLAB (R2021a or later recommended)

---

## 📌 Notes

This project was developed as part of a Navigation, Guidance & Control (GNC) evaluation task involving robust nonlinear re-entry guidance under uncertainty.

---

## 👤 Author

[Your Name]
