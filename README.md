# Separating thermal drift from heat loss by flow reversal

[Read the research manuscript (PDF)](paper/article.pdf) · [Numerical supplement](numerics/README-supplement.txt)

Can the attenuation of a thermal wave distinguish moving heat from heat lost to the walls? In a uniform channel with unknown heater and receiver gains, different combinations of transport speed and loss rate give identical one-way records. Equal-magnitude flow reversal supplies a second observation that separates these effects if the relative receiver response stays stable.

This repository contains a research draft, its editable LaTeX source, and simulations. No physical instrument measurements are reported. OpenAI Codex assisted with the analysis, implementation and writing, as disclosed in the manuscript.

## Reproduce the calculations

From the project directory, using Python 3.12:

    python -m pip install -r requirements.txt
    python numerics/verify.py
    python numerics/validation.py
    python numerics/study.py
    python numerics/regimes.py
    python numerics/artwork.py

The scripts write tables and figures beside the scripts. The paper directory contains the corresponding figures used to compile the manuscript. To rebuild after changing the figures, copy the regenerated layout.pdf and maps.pdf from numerics to paper, then run tectonic article.tex inside paper.

The baseline comparison uses the same 300 simulated observation sets for each estimator. The amplitude inversion recovered all three parameters within 10% in 16 trials. The complex fit with unknown receiver gains did so in 266 trials. See the paper for the noise model, the effect of gain drift, and conditions where accuracy declines.
