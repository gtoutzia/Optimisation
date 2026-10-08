# Optimisation Techniques

Implementations of classical optimisation algorithms from scratch, with a study of how their parameters affect convergence.
Course projects for *Optimisation Techniques*, Aristotle University of Thessaloniki.

## Project 1 – One-dimensional minimisation
Bisection, bisection with derivatives, golden-section and Fibonacci search (`Optimization1–4.m`).
Compares the number of function evaluations needed for different tolerances and search-interval lengths.

## Project 2 – Unconstrained multivariable minimisation
Minimisation of f(x,y) = x³·e^(−x²−y⁴) using:
- Steepest descent
- Newton's method
- Levenberg–Marquardt

Each method is tested with three step-size rules: constant step, step chosen by line search, and the Armijo rule.

## Project 3 – Constrained minimisation
Steepest descent with projection onto a feasible box, tested with different step sizes and starting points to show convergence and oscillation.

## Tools
MATLAB, Symbolic Math Toolbox (gradients and Hessians)

## How to run
Open MATLAB in the project folder and run the scripts. Each one prints the minimum found and plots convergence.

## Reports
Full reports (in Greek) are included as PDFs in each project folder.
