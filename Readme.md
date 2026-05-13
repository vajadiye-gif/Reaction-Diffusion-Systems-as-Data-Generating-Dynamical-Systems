# Reaction–Diffusion Systems as Data-Generating Dynamical Systems

> **Three reaction–diffusion models — Gray–Scott, Multi-fate, and EMT — explored as nonlinear physical systems and high-dimensional spatiotemporal data generators. Interactive Colab simulations with ipywidgets, parameter studies, and emergent pattern analysis.**

---

## Table of Contents

1. [Overview](#overview)
2. [Repository Structure](#repository-structure)
3. [Models](#models)
   - [Gray–Scott](#1-grayscott-model)
   - [Multi-fate](#2-multi-fate-model)
   - [EMT (Epithelial–Mesenchymal Transition)](#3-emt-model)
4. [Numerical Methods](#numerical-methods)
5. [Interactive Simulations](#interactive-simulations)
6. [Parameter Reference](#parameter-reference)
7. [Key Observations](#key-observations)
8. [Installation & Usage](#installation--usage)
9. [Notebooks Guide](#notebooks-guide)
10. [Future Work](#future-work)
11. [References](#references)

---

## Overview

Reaction–diffusion (RD) systems are PDEs of the form

$$
\frac{\partial \mathbf{u}}{\partial t} = \mathbf{D} \nabla^2 \mathbf{u} + \mathbf{f}(\mathbf{u})
$$

where **D** is a diagonal diffusion matrix and **f(u)** encodes the local nonlinear reaction kinetics. Despite their deceptively simple structure, RD systems generate a remarkable diversity of spatiotemporal behavior: stationary Turing patterns, travelling waves, spiral waves, and multi-stable domain coexistence.

This project treats these systems from two complementary perspectives:

**As physical models** — the patterns and steady states are interpreted in their biological and chemical context (chemical oscillators, cell fate decisions, cancer invasion).

**As data-generating dynamical systems** — each simulation produces a high-dimensional time series of spatial fields, which can be analysed with tools from nonlinear dynamics, time-series analysis, and data science (PCA of field snapshots, correlation structure, attractor reconstruction).

The three models span a progression of biological complexity:

| Model | Species | Phenomenon | Boundary Conditions |
|---|---|---|---|
| Gray–Scott | 2 (u, v) | Turing instability, self-replicating spots, labyrinths | Periodic |
| Multi-fate | 2 (A, B) | Bistable cell fate, domain coarsening, fate boundaries | Periodic / Neumann |
| EMT | 7 (miR-200, mZEB, ZEB, SNAIL, mSNAIL, miR-34, I) | Cancer invasion front, tristable gene circuit, spatial phenotype selection | Neumann |

---

## Repository Structure

```
Reaction-Diffusion-Systems-as-Data-Generating-Dynamical-Systems/
│
├── notebooks/                  # Interactive Jupyter notebooks (run in Colab)
│   ├── Gray_Scott.ipynb
│   ├── MultiFate.ipynb
│   └── EMT.ipynb
│
├── videos/                     # Simulation GIFs and MP4s
│   ├── Gray_Scott/
│   │   └── Labirynths (F = 0.0545, k=0.062).gif
│   ├── MultiFate/
│   │   └── Radial Gradient.gif
│   └── EMT/
│       └── SNAIL.gif
│
├── Images/                     # Static output figures
│
├── Readme.md
└── requirements.txt
```

> **Note on notebooks:** Due to interactive widgets and long-running simulations, notebooks are also exported as PDFs for static viewing on GitHub. The original interactive versions are intended for Google Colab.

---

## Models

### 1. Gray–Scott Model

The Gray–Scott model is the canonical activator-substrate reaction–diffusion system, originally proposed to describe autocatalytic chemical reactions. It is now widely used as a testbed for Turing instability and pattern formation theory.

#### Reaction Scheme

The model captures two processes: a self-replicating reaction where substrate $u$ is converted to product $v$ in the presence of $v$ itself (autocatalysis), and removal of both species at controlled rates.

$$
U + 2V \rightarrow 3V \qquad V \rightarrow P
$$

#### PDE System

$$
\begin{aligned}
\frac{\partial u}{\partial t} &= D_u \nabla^2 u - uv^2 + F(1-u) \\
\frac{\partial v}{\partial t} &= D_v \nabla^2 v + uv^2 - (F+k)v
\end{aligned}
$$

#### Parameters

| Symbol | Meaning | Typical Range |
|---|---|---|
| $D_u$ | Diffusion coefficient of substrate $u$ | 0.16 – 0.20 |
| $D_v$ | Diffusion coefficient of product $v$ | 0.04 – 0.08 |
| $F$ | Feed rate — how fast $u$ is replenished | 0.01 – 0.08 |
| $k$ | Kill rate — how fast $v$ is removed | 0.04 – 0.07 |

The ratio $D_u / D_v > 1$ is the necessary condition for Turing instability: the inhibitor ($v$) must diffuse more slowly than the activator ($u$).

#### Pattern Zoo

The $(F, k)$ parameter plane is a phase diagram of pattern morphologies. The labyrinth simulation shown here uses $F = 0.0545$, $k = 0.062$:

![Gray–Scott Labyrinth Pattern](https://github.com/vajadiye-gif/Reaction-Diffusion-Systems-as-Data-Generating-Dynamical-Systems/raw/main/videos/Gray_Scott/Labirynths%20(F%20%3D%200.0545%2C%20k%3D0.062).gif)

*Labyrinthine stripes at F = 0.0545, k = 0.062. The system evolves from a randomly perturbed homogeneous state; the Turing instability amplifies spatial fluctuations into stable, stationary stripe domains.*

> See `videos/Gray_Scott/` and `Images/` for the full parameter survey.

**Key observation:** Increasing $F$ at fixed $k$ accelerates the coarsening of pattern domains and shifts the morphology from labyrinthine stripes toward isolated spots. Near the boundary between stable patterns and the trivial steady state, the system exhibits transient irregular dynamics before settling.

---

### 2. Multi-fate Model

The Multi-fate model captures bistable cell fate decisions in a developing tissue. Species $A$ and $B$ represent two mutually repressive transcription factors (or gene programs), each positively autoregulated through dimerization. The spatial coupling via diffusion allows different regions of the tissue to settle into different fate domains, producing a patterned tissue from a nearly homogeneous initial condition.

#### PDE System

$$
\begin{aligned}
\frac{\partial A}{\partial t} &= D_u \nabla^2 A + f_A(A, B) \\
\frac{\partial B}{\partial t} &= D_v \nabla^2 B + f_B(A, B)
\end{aligned}
$$

#### Reaction Terms

The reaction functions encode autoactivation through dimerization and basal production:

$$
\begin{aligned}
f_A(A,B) &= \alpha + \frac{\beta \, t_A}{1 + t_A} - A \\
f_B(A,B) &= \alpha + \frac{\beta \, t_B}{1 + t_B} - B
\end{aligned}
$$

where $\alpha$ is the basal production rate and $\beta$ is the maximum activation strength.

#### Dimerization and Hill Activation

The cooperative activation terms $t_A$, $t_B$ arise from the equilibrium dimer concentrations $A_2$, $B_2$:

$$
t_A = (A_2)^n, \qquad t_B = (B_2)^n
$$

$$
\begin{aligned}
A_2 &= \frac{2A^2}{K + 4(A + B) + \sqrt{K^2 + 8K_d(A + B)}} \
B_2 &= \frac{2B^2}{K + 4(A + B) + \sqrt{K^2 + 8K_d(A + B)}}
\end{aligned}
$$

The shared denominator reflects **competitive dimerization**: $A$ and $B$ monomers compete for the same dimerization pool, coupling the two species nonlinearly even before diffusion.

| Symbol | Meaning |
|---|---|
| $K$ | Total monomer pool size |
| $K_d$ | Dimerization dissociation constant |
| $n$ | Hill coefficient for cooperative activation |
| $\alpha$ | Basal production rate |
| $\beta$ | Maximum activation gain |

#### Sample Output

![MultiFate Radial Pattern](https://github.com/vajadiye-gif/Reaction-Diffusion-Systems-as-Data-Generating-Dynamical-Systems/raw/main/videos/MultiFate/Radial%20Gradient.gif)

*Radial gradient initial condition: a central seed of fate-A state surrounded by fate-B background. The interface propagates radially as the system resolves the spatial competition between the two stable attractors.*

> See `videos/MultiFate/` and `Images/` for additional initial conditions and parameter regimes.

**Key observation:** Different diffusion constant ratios $D_u / D_v$ produce qualitatively distinct spatial outcomes — when diffusion is symmetric ($D_u = D_v$), fate domains coarsen by curvature-driven interface motion (Allen–Cahn dynamics). When the ratio is large, the faster-diffusing species can outrun the other and invade larger territories, breaking the symmetry of domain sizes.

---

### 3. EMT Model

The EMT model is the spatially-extended version of the SNAIL/ZEB/miR-200/miR-34 gene regulatory circuit, which governs the Epithelial–Mesenchymal Transition — a central mechanism in cancer metastasis and embryonic development. Diffusion here models the intercellular spread of signaling molecules across a 2D tissue.

#### Circuit Biology

The circuit is built from two mutually repressing double-negative feedback loops:

```
miR-200 ──┤ ZEB   ──┤ miR-200    (bistable switch)
miR-34  ──┤ SNAIL ──┤ miR-34     (bistable switch)
SNAIL   ──→ ZEB                  (positive feed-forward coupling)
```

Together these create a **tristable** system at the single-cell level, with three coexisting attractors: Epithelial (E), Hybrid E/M, and Mesenchymal (M).

#### PDE System (7 species)

$$
\begin{aligned}
\frac{\partial [\text{miR200}]}{\partial t} &= G_\text{MIR200} \cdot h_\text{mir200}^\text{zeb} \cdot h_\text{mir200}^\text{sna} - [\text{miR200}]^1 \cdot \text{mir200}_0 - K_\text{MIR200}\,[\text{miR200}] + D_\text{mir200}\,\nabla^2[\text{miR200}] \\[4pt]
\frac{\partial [m_\text{ZEB}]}{\partial t} &= G_\text{MZEB} \cdot h_\text{mzeb}^\text{zeb} \cdot h_\text{mzeb}^\text{sna} - [\text{miR200}]^1 \cdot \text{mir200}_1 - K_\text{MZEB}\,[m_\text{ZEB}] + D_\text{mzeb}\,\nabla^2[m_\text{ZEB}] \\[4pt]
\frac{\partial [\text{ZEB}]}{\partial t} &= G_\text{ZEB}\,[\text{miR200}]^1 \cdot \text{mir200}_2 - K_\text{ZEB}\,[\text{ZEB}] + D_\text{ZEB}\,\nabla^2[\text{ZEB}] \\[4pt]
\frac{\partial [\text{SNAIL}]}{\partial t} &= G_\text{SNAIL}\,[\text{miR34}]^4 \cdot \text{mir34}_2 - K_\text{SNAIL}\,[\text{SNAIL}] + D_\text{SNAIL}\,\nabla^2[\text{SNAIL}] \\[4pt]
\frac{\partial [M_\text{SNAIL}]}{\partial t} &= G_\text{MSNAIL} \cdot h_\text{msnail} \cdot h_\text{msna}^\text{sna} - [\text{miR34}]^4 \cdot \text{mir34}_1 - K_\text{MSNAIL}\,[M_\text{SNAIL}] + D_\text{MSNAIL}\,\nabla^2[M_\text{SNAIL}] \\[4pt]
\frac{\partial [\text{miR34}]}{\partial t} &= G_\text{MIR34} \cdot h_\text{mir34}^\text{zeb} \cdot h_\text{mir34}^\text{sna} - [\text{miR34}]^4 \cdot \text{mir34}_0 - K_\text{MIR34}\,[\text{miR34}] + D_\text{mir34}\,\nabla^2[\text{miR34}] \\[4pt]
\frac{\partial [I]}{\partial t} &= 0
\end{aligned}
$$

where $I$ is the external TGF-$\beta$ proxy signal, spatially pinned (no diffusion, no dynamics).

The inhibitory Hill functions $h$ encode transcriptional repression:

$$
h_x^y = \frac{1}{1 + \left(\frac{[y]}{\theta_{xy}}\right)^{n_{xy}}}
$$

and the miRNA–mRNA binding terms (e.g. $\text{mir200}_0$, $\text{mir200}_1$, $\text{mir200}_2$) arise from a combinatorial binding model summed over the occupation states of the multiple miRNA binding sites on each mRNA.

#### Phenotype Classification

The cellular phenotype $\mathcal{P}$ at grid point $(i,j)$ is classified from the steady-state miR-200/ZEB ratio:

$$
\mathcal{P}_{i,j} =
\begin{cases}
\text{Epithelial (E)} & \text{if } [\text{miR200}] \gg [\text{ZEB}] \\
\text{Hybrid (E/M)} & \text{if } [\text{miR200}] \approx [\text{ZEB}] \\
\text{Mesenchymal (M)} & \text{if } [\text{miR200}] \ll [\text{ZEB}]
\end{cases}
$$

#### Phenotypic Molecular Signatures

$$
\begin{array}{|l|c|c|c|}
\hline
\textbf{Phenotype} & [\text{miR200}] & [\text{ZEB}] & \text{Biological State} \\
\hline
\color{blue}{\text{Epithelial (E)}} & \uparrow\,\text{High} & \downarrow\,\text{Low} & \text{Adherent / Stationary} \\
\hline
\color{green}{\text{Hybrid (E/M)}} & \sim\,\text{Mid} & \sim\,\text{Mid} & \text{Collective Migration} \\
\hline
\color{red}{\text{Mesenchymal (M)}} & \downarrow\,\text{Low} & \uparrow\,\text{High} & \text{Invasive / Solitary} \\
\hline
\end{array}
$$

#### Sample Output

![EMT SNAIL Pattern](https://github.com/vajadiye-gif/Reaction-Diffusion-Systems-as-Data-Generating-Dynamical-Systems/raw/main/videos/EMT/SNAIL.gif)

*SNAIL concentration field evolving on the 2D tissue grid. A mesenchymal-state seed expands into the epithelial background as the invasion front propagates. The spatial sharpness of the E/M interface depends on the relative diffusion rates.*

> See `videos/EMT/` and `Images/` for ZEB, miR-200, and miR-34 field animations.

**Key observation:** When the external signal $I$ is in the bistable regime (between the two saddle-node bifurcation points), the M-state patch expands while the E-state background remains stable. The Hybrid E/M branch, which exists in the single-cell ODE model, is destroyed by spatial diffusion — the tissue-level system is bistable even when the cell-level system is tristable.

---

## Numerical Methods

All three models are integrated on a 2D spatial grid using the same core scheme:

**Spatial discretisation:** 5-point finite difference stencil for the Laplacian on a uniform Cartesian grid with spacing $\Delta x = 1$:

$$
\nabla^2 u_{i,j} \approx \frac{u_{i+1,j} + u_{i-1,j} + u_{i,j+1} + u_{i,j-1} - 4u_{i,j}}{\Delta x^2}
$$

**Time integration:** Explicit Euler for diffusion; the reaction step is integrated with either explicit Euler or RK4 depending on the model stiffness.

**Boundary conditions:**

| Model | BC Type | Implementation |
|---|---|---|
| Gray–Scott | Periodic | Array wrap (`np.roll`) |
| Multi-fate | Periodic / Neumann | Ghost-node method for Neumann |
| EMT | Neumann (zero-flux) | Ghost-node: boundary cell's missing neighbour replaced by the nearest interior cell, enforcing $\partial u/\partial n = 0$ |

**Stability (CFL condition for explicit diffusion):**

$$
\frac{D_\text{max} \cdot \Delta t}{\Delta x^2} \leq \frac{1}{2}
$$

Timesteps are chosen to satisfy this bound with a safety margin.

---

## Interactive Simulations

Notebooks are built for Google Colab with `ipywidgets` for real-time interactivity:

- **ToggleButtons** — switch the displayed species field (e.g., $u$/$v$ for Gray–Scott; SNAIL/ZEB/miR-200/miR-34 for EMT)
- **Play/IntSlider** — scrub through saved simulation frames with animated playback
- **Parameter sliders** — adjust $F$, $k$, diffusion coefficients, and signal level $I$ and rerun simulations interactively

```python
# Colab setup — run once per session
from google.colab import output
output.enable_custom_widget_manager()  # required for Play button in Colab
```

---

## Parameter Reference

### Gray–Scott

| Parameter | Labyrinth | Spots | Traveling Waves |
|---|---|---|---|
| $F$ | 0.0545 | 0.035 | 0.022 |
| $k$ | 0.062 | 0.060 | 0.051 |
| $D_u$ | 0.16 | 0.16 | 0.16 |
| $D_v$ | 0.08 | 0.08 | 0.08 |

### Multi-fate

| Parameter | Symbol | Role |
|---|---|---|
| Basal production | $\alpha$ | Sets the baseline activity of both fates |
| Max activation | $\beta$ | Controls how strongly dimerization drives autoactivation |
| Hill coefficient | $n$ | Sharpness of the bistable switch |
| Dimerization constant | $K_d$ | Controls the coupling between A and B through the shared dimer pool |

### EMT

See the full kinetic parameter table in the companion [EMT-Circuit-Analysis](https://github.com/vajadiye-gif/EMT-Circuit-Analysis) repository.

---

## Key Observations

**Gray–Scott:** The $(F, k)$ parameter plane is a morphology phase diagram. Labyrinths, spots, worms, and homogeneous states occupy distinct regions; parameter boundaries exhibit hysteresis. Increasing $F$ at fixed $k$ accelerates spatial coarsening.

**Multi-fate:** Asymmetric diffusion ($D_u \neq D_v$) breaks the symmetry of fate domain sizes, with the faster-diffusing species occupying a larger territory. Symmetric diffusion produces curvature-driven interface dynamics (Allen–Cahn coarsening). The fate boundary sharpness is controlled by $K_d$.

**EMT:** The tissue-level bistability is a consequence of diffusion suppressing the Hybrid E/M branch that exists in the single-cell model. The M-state invasion front speed depends on the signal level $I$ — at the lower saddle-node point, front propagation stalls. Varying reaction rates changes cluster morphology from sharp invasion fronts to diffuse mixed domains.

---

## Installation & Usage

```bash
# Clone the repository
git clone https://github.com/vajadiye-gif/Reaction-Diffusion-Systems-as-Data-Generating-Dynamical-Systems.git
cd Reaction-Diffusion-Systems-as-Data-Generating-Dynamical-Systems

# Install dependencies
pip install -r requirements.txt
```

### Google Colab (recommended for interactivity)

```python
# Mount Drive (optional — for saving GIFs/arrays)
from google.colab import drive
drive.mount('/content/drive')

# Clone if not already present
import os
if not os.path.isdir('Reaction-Diffusion-Systems-as-Data-Generating-Dynamical-Systems'):
    !git clone https://github.com/vajadiye-gif/Reaction-Diffusion-Systems-as-Data-Generating-Dynamical-Systems.git

# Enable widgets
from google.colab import output
output.enable_custom_widget_manager()
```

Then open any notebook in `notebooks/`:

1. Select the species to visualize using the **ToggleButtons** widget
2. Run all cells — the simulation will execute and stream frames to the interactive heatmap
3. Use the **Play/Slider** widget to animate the saved frames
4. Save outputs as `.npy` arrays, `.png` images, or `.gif` videos as needed

---

## Notebooks Guide

### `Gray_Scott.ipynb`
Simulates the Gray–Scott system on a 2D periodic grid. Includes parameter sweep widgets and animated visualisation of $u$ and $v$ fields. Demonstrates Turing instability from random initial perturbations.

### `MultiFate.ipynb`
Simulates the Multi-fate bistable circuit with several initial conditions (uniform, radial gradient, random). Visualises species $A$ and $B$ and the emergent fate domain structure. Explores how diffusion ratio alters domain size distributions.

### `EMT.ipynb`
Simulates the 7-species EMT circuit on a 2D grid with Neumann BCs. Initialises an M-state patch in an E-state background. Animates all 6 dynamic species fields (miR-200, mZEB, ZEB, SNAIL, mSNAIL, miR-34). Demonstrates spatial bistability and M-state invasion.

---

## Future Work

- **3D pattern formation:** Extend the Laplacian to a 3D grid and initialise a spherical M-core (tumour spheroid geometry)
- **Turing instability analysis:** Linear stability analysis of the homogeneous steady state to analytically predict the onset wavenumber $k^*$ and compare with observed pattern wavelengths
- **Stochastic extension:** Add Langevin noise to the reaction terms to model intrinsic molecular fluctuations near bifurcation points
- **Parameter sensitivity / bifurcation sweeps:** Automated sweeps over $(F, k)$ for Gray–Scott and $(D_u/D_v)$ for Multi-fate to map full morphology phase diagrams
- **Data science perspective:** PCA of field snapshots, attractor reconstruction from time series, and anomaly detection in parameter-driven pattern transitions

---

## References

1. Gray, P. & Scott, S. K. (1985). Sustained oscillations and other exotic patterns of behavior in isothermal reactions. *J. Phys. Chem.*, 89(1), 22–32.
2. Maire, T. & Youk, H. (2015). Molecular-level tuning of cellular autonomy controls the collective behaviors of cell populations. *Cell Systems*, 1(5), 349–360. *(Multi-fate model)*
3. Lu, M., Jolly, M. K., Levine, H., Onuchic, J. N., & Ben-Jacob, E. (2013). MicroRNA-based regulation of epithelial–hybrid–mesenchymal fate determination. *PNAS*, 110(45), 18144–18149. *(EMT circuit)*
4. Cross, M. C. & Hohenberg, P. C. (1993). Pattern formation outside of equilibrium. *Reviews of Modern Physics*, 65(3), 851.

---

**Supervisor:** Ushasi Roy, IISER Pune  
**Course:** Semester project in computational biophysics, 2025–2026  
**Author:** Ved Amar Jadiye (BSMS Physics, IISER Pune)

Built with **NumPy · Matplotlib · ipywidgets · Numba · Google Colab**
