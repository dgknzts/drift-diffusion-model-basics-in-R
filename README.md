# Diffusion Decision Model (DDM) Basics in R - From Scratch

This repository provides a resource for understanding and implementing the fundamental principles of the Diffusion Decision Model (DDM) from scratch using the R programming language. The goal is to present the underlying logic and simulation of the model in the simplest and most accessible way.

This work is heavily inspired by and aims to replicate concepts from "An Introduction to the Diffusion Model of Decision-Making" by Philip L. Smith and Roger Ratcliff.

## Repository Structure

*   **R Functions (`R/`):**
    *   Simulators for simple Random Walks and Diffusion Models
    *   Functions to demonstrate the effects of DDM parameters
    *   Helper functions for visualizing results

*   **Documentation (`vignettes/`):**
    *   Introduction to random walks and their simulation
    *   Simulating the diffusion model and its logic
    *   Understanding DDM parameters (drift rate, threshold, starting point, etc.)
    *   Exploring parameter variability in DDM


*   **Shiny App (`local-app/`):**
    *   Interactive app for exploring DDM parameters and visualizing results
    *   Published online at: https://dgknzts.shinyapps.io/DDM-Basics/
    *   Basic demonstration version of the DDM Shiny app for easy local execution
    *   Simplified version for educational purposes

*   **Documentation Site (`docs/`):**
    *   Generated pkgdown site with tutorials and API reference

## Getting Started

1. Clone this repository: `git clone https://github.com/dgknzts/drift-diffusion-model-basics-in-R.git`
2. Open `DDM_Basics_R.Rproj` in RStudio
3. Explore the vignettes to learn about DDM:
   * Visit the documentation site at https://dgknzts.github.io/drift-diffusion-model-basics-in-R/
   * Or open the R Markdown files in the `vignettes/` folder

## Interactive DDM Explorer

To run the interactive Shiny app locally:

```r
# In R console
shiny::runApp("local-app")
```

Or visit the hosted version at: https://dgknzts.shinyapps.io/DDM-Basics/

## Key DDM Parameters Implemented

*   `v`: Drift rate (mean rate of evidence accumulation)
*   `a`: Threshold separation (boundary separation)
*   `z`: Starting point of evidence accumulation (bias)
*   `ter` (or `t0`): Non-decision time (encoding and response execution time)
*   `s`: Within-trial noise (standard deviation of the normally distributed increments)
*   Parameter variability: `sv`, `sz`, `st0` (across-trial variability in drift, starting point, and non-decision time)
    
## References

*   Smith, P. L., & Ratcliff, R. (2024). An Introduction to the Diffusion Model of Decision-Making. *In B. U. Forstmann & B. M. Turner (Eds.), An Introduction to Model-Based Cognitive Neuroscience.*
*   Ratcliff, R. (1978). A theory of memory retrieval. *Psychological Review, 85*(2), 59–108.

## License

This project is licensed under the MIT License - see the `LICENSE` file for details.
