Numerical supplement | Harmonic flow reversal

All numbers are synthetic; lengths in mm, time in s. No apparatus or measurement files are represented. Python 3.12 and NumPy 2.x were used. Plotting additionally requires ReportLab. Run from an extracted directory. No downloads are performed by the scripts.

python verify.py
  One-dimensional harmonic finite-difference solver and two-frequency inversion. Recreates results.csv, the clean-grid and two 3000-realization noise checks.
python validation.py
  Checks raw complex responses against the finite-domain Dirichlet Green function. Recreates grid.csv and prints results for alternative optimizer starts.
python study.py
  Produces benchmark.csv from paired simulations. “complex” estimates the receiver ratio at each frequency. “calibrated” supplies its true nominal value. When the receiver drifts on reversal, even this nominal calibration is misspecified.
python regimes.py
  Produces regimes.csv for the two-dimensional operating map, with 200 realizations per cell.
python artwork.py
  Recreates the two vector illustrations from the geometry and the map CSV. Requires ReportLab.

Data dictionary
noise = complex rms of the synthetic circular Gaussian error in response units.
minSNR = smallest clean response amplitude / noise.
physical = number out of 300 with finite estimates, positive D and nonnegative loss.
within10 = number out of all 300 with physical estimates and each relative error <0.10.
medD, medV, medLoss = median among physical estimates.
MdARE_D, MdARE_v, MdARE_loss = median absolute relative errors (fractions) among physical estimates.
Pe = vL/D. loss_number = lambda L^2/D. Map entries are counts out of 200.

The benchmark and map use common Gaussian draws across settings, with different seeds for the two experiments. A repetition starts in the same deterministic noiseless forward state; only the observation noise is resampled. The random seeds appear in each driver. The source files contain the full damping schedule and parameter bounds. The underlying simulated observations can be regenerated from these inputs.
