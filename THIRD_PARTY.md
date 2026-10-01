# Dependencies and attribution

This release contains research code and simulated results. It does not redistribute MATLAB, CasADi runtime binaries, Python packages, third-party paper PDFs or publisher templates. The `.casadi` files are serialized project model functions and require an independently installed compatible CasADi runtime.

Install dependencies from their official distributions. Their own terms continue to apply independently of this project's MIT license:

- [CasADi](https://web.casadi.org/): nonlinear modeling, automatic differentiation, integrators and optimization interfaces. See its [documentation and citation guidance](https://web.casadi.org/docs/).
- [MATLAB](https://www.mathworks.com/products/matlab.html): separately licensed runtime for the MATLAB model and ARC runners.
- [NumPy](https://numpy.org/), [SciPy](https://scipy.org/), [Matplotlib](https://matplotlib.org/), [SymPy](https://www.sympy.org/) and [PyYAML](https://pyyaml.org/): install separately using `requirements.txt`.

CasADi distributions can include further numerical solver components. Their notices belong to those distributions; this project does not relicense them. Scientific references for the reactor equations and control analysis belong to the associated manuscript. Implementing or citing an equation is distinct from claiming authorship of a third-party library.
