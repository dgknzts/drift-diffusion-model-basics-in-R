# Diffusion Decision Model (DDM) Basics in R - From Scratch

This repository provides a resource for understanding and implementing the fundamental principles of the Diffusion Decision Model (DDM) from scratch using the R programming language. The goal is to present the underlying logic and simulation of the model in the simplest and most accessible way.

This work is heavily inspired by and aims to replicate concepts from "An Introduction to the Diffusion Model of Decision-Making" by Philip L. Smith and Roger Ratcliff.

## Contents

*   **R Functions (`R/`):**
    *   Simulators for simple Random Walks.
    *   Simulators for the basic Diffusion Model.
    *   Functions to demonstrate the effects of DDM parameters.
    *   Helper functions for visualizing results.
*   **Vignettes/Tutorials (`vignettes/`):**
    *   Introduction to random walks and their simulation.
    *   Simulating the diffusion model and its logic.
    *   Understanding DDM parameters (drift rate, threshold, starting point, etc.) and their impact.
    *   How to plot DDM outputs (RT distributions, choice probabilities).
*   **Example Scripts (`scripts/`):**
    *   Running simple simulations using the provided functions.
    *   Generating example plots.

## Getting Started

1.  Clone this repository: `git clone https://github.com/yourusername/DDM_Basics_R.git`
2.  Open `DDM_Basics_R.Rproj` in RStudio.
3.  Start by exploring the R Markdown files in the `vignettes/` folder.
4.  Run the example scripts in the `scripts/` folder to see how the functions are used.

## Key DDM Parameters Implemented

*   `v`: Drift rate (mean rate of evidence accumulation)
*   `a`: Threshold separation (boundary separation)
*   `z`: Starting point of evidence accumulation (bias)
*   `ter` (or `t0`): Non-decision time (encoding and response execution time)
*   `s`: Within-trial noise (standard deviation of the normally distributed increments, often fixed to 0.1)
    *   *(Parameter variability like `sv`, `sz`, `st` can be added as an extension)*

## Contributing

Contributions and feedback are welcome! Please open an issue or submit a pull request.

## References

*   Smith, P. L., & Ratcliff, R. (YEAR). An Introduction to the Diffusion Model of Decision-Making. *In B. U. Forstmann & B. M. Turner (Eds.), An Introduction to Model-Based Cognitive Neuroscience.* (This is the primary chapter this repo aims to illustrate).
*   Ratcliff, R. (1978). A theory of memory retrieval. *Psychological Review, 85*(2), 59–108.

## License

This project is licensed under the MIT License - see the `LICENSE` file for details.