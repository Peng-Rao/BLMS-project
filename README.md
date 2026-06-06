# Bayesian Learning & Monte Carlo Simulation — Final Project 2026

## Yacht Hydrodynamic Resistance

A Bayesian study of the residuary resistance of sailing yachts as a function of
hull geometry and speed (Froude number).

## Dataset

**Source.** UCI Machine Learning Repository,
<https://archive.ics.uci.edu/dataset/243/yacht+hydrodynamics>.
Original record: Gerritsma, Onnink and Versluis, *Geometry, resistance and
stability of the Delft Systematic Yacht Hull Series*, 1981.

**Data.** The dataset describes 22 different hull forms (sailing yachts), each
tested at multiple Froude numbers in a towing tank, producing 308 hull-form ×
Froude-number rows (`data/yacht_hydro.csv`).

**Variables.**

| Variable       | Description                                                        |
| -------------- | ----------------------------------------------------------------- |
| `LCB`          | Longitudinal position of the center of buoyancy (% length)        |
| `prismatic`    | Prismatic coefficient (dimensionless)                             |
| `L_disp`       | Length–displacement ratio                                         |
| `beam_draught` | Beam–draught ratio                                               |
| `L_beam`       | Length–beam ratio                                                |
| `Froude`       | Froude number (dimensionless, *Fn = v/√(gL)*), in [0.125, 0.45]   |
| `resistance`   | Residuary resistance per unit weight of displacement (**target**) |

The first five covariates are geometric form factors that are constant within a
given hull. `Froude` is a dimensionless measure of speed and is the dominant
predictor: for sailing yachts the residuary resistance scales approximately as a
power law of the Froude number, with an exponent close to 4 in the bulk regime.

## Task

The physics of the problem says that for sailing yachts the residuary resistance
follows a power law of the Froude number with exponent close to 4 in the bulk
regime. The project explores this with a Bayesian analysis:

1. **Log–log regression.** Fit

   $$\log(\mathrm{resistance}) = \beta_0 + \beta_F\,\log(\mathrm{Froude}) + \sum_k \beta_k\,\mathrm{hull}_k + \varepsilon\,.$$

   Verify that $\hat\beta_F \approx 4$.

2. **BAS / BMA.** Use `bas.lm` to enumerate the $2^5$ submodels over the five
   hull-geometry covariates and identify which of them (LCB, prismatic, L/disp,
   B/draught, L/beam) carry signal once `Froude` is in the model.

3. **Correct for heteroscedasticity.** Train the nonlinear model

   $$\mathrm{resistance} = k\,\mathrm{Froude}^{\gamma_F}\,\exp\left\lbrace\sum_k \beta_k\,\mathrm{hull}_k\right\rbrace + \varepsilon$$

   where only the form factors that appear relevant in the model-selection phase
   are kept. Compare the results with the log–log regression and discuss the
   impact of the heteroscedasticity correction on the inference on $\gamma_F$ and
   on the hull-geometry covariates.

4. **Posterior predictive.** Plot the posterior predictive curves and compare the
   curves with and without the treatment of the heteroscedasticity, and compare
   the curves with the data.

## Repository layout

| Path                     | Contents                                              |
| ------------------------ | ----------------------------------------------------- |
| `data/yacht_hydro.csv`   | The yacht hydrodynamics dataset (UCI #243)            |
| `yacht_resistance.Rmd`   | Analysis notebook (R Markdown)                        |
| `report.tex` / `report.pdf` | Project report (LaTeX source and compiled PDF)     |
| `figures/`               | Figures used in the report                            |
| `Project2026_INFO.pdf`   | Official final-project instructions                   |
