# Quantum Mechanics Projects - Eigenstates, Time-Dependent Solutions & Numerical Solvers

The dynamics of a one-dimensional quantum system are governed by the time-dependent Schrödinger equation:

$$ i\hbar\frac{\partial \psi}{\partial t} = -\frac{\hbar^2}{2m}\frac{\partial^2 \psi}{\partial x^2} + V\psi $$

The wave function $\psi$ is a function of both position $x$ and time $t$, and is the fundamental description of the realm of the very small. Imagine we are following the motion of a single particle in one dimension. This wave function represents a probability of measuring the particle at a position $x$ at a time $t$. Quantum mechanics tells us that (contrary to our familiar classical reasoning) this probability is not a limitation of our knowledge of the system, but a reflection of an unavoidable uncertainty about the position and time of events in the realm of the very small.

To solve this equation numerically, we'll use the split-step Fourier method.

## The Split-Step Fourier Method

A common way to numerically solve certain differential equations is through the use of the Fourier transform. We'll use a Fourier convention of the following form:

$$ \tilde{\psi}(k, t) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{\infty} \psi(x, t) e^{-ikx} dx $$

Under this convention, the associated inverse Fourier Transform is given by:

$$ \psi(x, t) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{\infty} \tilde{\psi}(k, t) e^{ikx} dk $$

Substituting this into the Schrödinger equation and simplifying gives the Fourier-space form of the Schrödinger equation:

$$ i\hbar\frac{\partial \tilde{\psi}}{\partial t} = -\frac{\hbar^2}{2m}k^2\tilde{\psi}+V\psi $$

The two versions of the Schrödinger equation contain an interesting symmetry: the time step in each case depends on a straightforward multiplication of the wave function $\psi$, as well as a more complicated term involving derivatives with respect to $x$ or $k$. The key observation is that while the equation is difficult to evaluate fully within one of the forms, each basis offers a straightforward calculation of one of the two contributions. This suggests an efficient strategy to numerically solve the Schrödinger Equation.

### Numerical Considerations

Solving this system numerically will require repeated computations of the Fourier transform of $\psi(x, t)$ and the inverse Fourier transform of $\tilde{\psi}(k, t)$. The most known algorithm for computing numerical Fourier transforms is the Fast Fourier Transform (FFT), available in NumPy and efficiently computing the following form of the discrete Fourier transform:

$$ F\tilde{\psi}= \sum_{n=0}^{N-1} f_n e^{-2\pi i n m/N} $$

where $f_n$ and $F\tilde{\psi}$ are the input and output vectors respectively. We need to know how these relate to the continuous Fourier transforms defined and used above. For instance, assume that the infinite integral is well-approximated by the finite integral from $a$ to $b$, so that we can write

$$ \tilde{\psi}(k, t) = \frac{1}{\sqrt{2\pi}}\int_{a}^{b} \psi(x, t) e^{-ikx} dx $$

This approximation ends up being equivalent to assuming that the potential $V(x) \rightarrow \infty$ at $x < a$ and $x > b$. We'll now approximate this integral as a Riemann sum of $N$ terms, and define $\Delta x = (b - a) / N$, and $x_n = a + n\Delta x$:

$$ \tilde{\psi}(k, t) \approx \frac{1}{\sqrt{2\pi}}\sum_{n=0}^{N-1} \psi(x_n, t) e^{-ikx_n} \Delta x $$

This is starting to look like the discrete Fourier transform! To bring it even closer, let's define $k_m = k_0 + m\Delta k$, where $\Delta k = 2\pi / (N\Delta x)$. Then our approximation becomes

$$ \tilde{\psi}(k_m, t) \approx \frac{1}{\sqrt{2\pi}}\sum_{n=0}^{N-1} \psi(x_n, t) e^{-ik_0x_n} e^{-2\pi imn/\mathcal{N}} \Delta x $$

(Note that just as we have limited the range of $x$ above, we have here limited the range of $k$ as well. This means that high-frequency components of the signal will be lost in our approximation. The Nyquist sampling theorem tells us that this is an unavoidable consequence of choosing discrete steps in space, and it can be shown that the spacing we chose above exactly satisfies the Nyquist limit if we choose $k_0 = -\pi / \Delta x$)

Plugging our expressions for $x_n$ and $k_m$ into the Fourier approximation and rearranging, we find the following:

$$ [\psi(x_n, t)e^{ik_0x_n}] \approx \sum_{m=0}^{\mathcal{N}-1}[\Delta x\sqrt{2\pi}\psi(x_m, t)e^{-2\pi imn/\mathcal{N}}]e^{ik_mx_n} $$

which matches the definition of the inverse DFT when we identify $\Delta x\sqrt{2\pi}$ as the weight factor for the coefficients. Therefore, the relationship between the continuous and discrete Fourier pairs can be expressed as follows:

$$ \psi(x, t) \rightleftarrows \psi(x_n, t) \\ \tilde{\psi}(k, t) \rightleftarrows \psi(x_n, t)\Rightarrow \tilde{\psi}(k_m, t) $$

under the specified conditions. This allows for a fast numerical evaluation of the Schrödinger equation.

### Putting It All Together: the Algorithm

Choose $a$, $b$, $N$, and $k_0$ such that they sufficiently represent the initial state of your wave function $\psi(x)$. (Caution: this might be the hardest part of the entire solution. If limits in $x$ or $k$ are chosen that don't fit your problem, the approximations we've made could ruin the accuracy of the calculation.) Once these are selected, then $\Delta x = (b - a) / N$ and $\Delta k = 2\pi / (N\Delta x)$. Define $x_n = a + n \Delta x$ and $k_m = k_0 + m \Delta k$.

Discretize the wave functions on this grid. Set $\psi_n(t) = \psi(x_n, t)$, $V_n = V(x_n)$, and $\tilde{\psi}_m = \tilde{\psi}(k_m, t)$.

To advance the system by a time step $\Delta t$, carry out the following:

* Compute a half-step in $x$: $\psi_n \gets \psi_n \exp[-i(\Delta t/2)(V_n/\hbar)]$
* Calculate $\tilde{\psi}_m$ from $\psi_n$ using the FFT.
* Perform a full-step in $k$: $\tilde{\psi}_m \gets \tilde{\psi}_m \exp[-i\hbar(k\cdot k)\Delta t/(2m)]$
* Calculate $\psi_n$ from $\tilde{\psi}_m$ using the inverse FFT.
* Carry out another half-step in $x$: $\psi_n \gets \psi_n \exp[-i(\Delta t/2)(V_n/\hbar)]$

Repeat step 3 until the desired time is achieved.

Please note that we have divided the $x$-space time-step into two half-steps; this leads to a more stable numerical solution compared to executing the whole time-step simultaneously. Some may recognize this as an implementation of the widely recognized leap-frog integration scheme.

As an illustrative example, I wrote Python code setting up a particle in a box with a potential barrier. The barrier was set high enough that a classical particle couldn't pass through, yet a quantum particle could still "tunnel," resulting in a non-zero chance of finding the particle on the opposite side of the barrier. This phenomenon, called quantum tunneling, underlies various technological applications, including electron microscopes, semiconductor diodes, and possibly future low-power transistors.

Below is the animated result (for a brief introduction to Python's animation capabilities, please check this article). The upper panel depicts the position-space wave function, while the lower panel displays the momentum-space wave function.

Even though the height of the potential barrier (indicated by the dashed line in the lower panel) exceeds the energy of the particle, due to quantum effects, some portion of the wavefunction manages to tunnel through the barrier and arrive on the other side.

## Table of Contents

- [Eigenstates of 1D Potential](https://github.com/tetraethylmethane/SchrodingerEquation/blob/main/Eigenstates%20of%201D%20Potential.ipynb)

To calculate eigenstates (energy levels) of a particle in a one-dimensional potential well or barrier using the numerical method of your choice. You can select from different types of potentials such as finite square well, infinite square well, harmonic oscillator, etc.

- [Schrödinger Equation Solution](https://github.com/tetraethylmethane/SchrodingerEquation/blob/main/Schr%C3%B6dinger%20Equation%20Solution.ipynb)

To find solutions to the stationary state Schrodinger equation for various potentials. The methods used include analytical solutions, shooting methods, and numerical integration techniques like the finite difference method.

- [Time-Dependent Schrödinger Equation](https://github.com/tetraethylmethane/SchrodingerEquation/blob/main/Time-Dependent%20Schrodinger%20Equation.ipynb)

To solve the time-dependent Schrodinger equation for different initial conditions and potentials. Techniques such as the split-operator method and Crank-Nicolson scheme are implemented to obtain the time evolution of wavefunctions.

- [General Numerical Solver for the 1D Time-Dependent Schrödinger's Equation Using the Split-Step Fourier Method (Schrödinger.py)](https://github.com/tetraethylmethane/SchrodingerEquation/blob/main/schrodinger.py)

Provides a general solver for the 1D time-dependent Schrodinger equation using the split-step Fourier method. It allows you to input custom potentials and initial conditions, making it versatile for various applications.

## Installation
To run these scripts and notebooks, you will need the following dependencies installed:

* Python 3.6+
* Jupyter Notebook/JupyterLab
* NumPy
* Matplotlib
* SciPy
* SymPy (optional, for symbolic computation)
* Plotly
* Torch
* scikit-image
* Science Plots
* Numba
* Sympy

You can install the necessary packages by running:
```bash
pip install numpy matplotlib scipy sympy numba SciencePlots scikit-image torch plotly
