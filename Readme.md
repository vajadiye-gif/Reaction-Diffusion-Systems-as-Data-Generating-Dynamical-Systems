# Reaction-Diffusion Project

This project simulates **reaction-diffusion systems** using Python. It includes three models: **Gray-Scott**, **Multi-fate**, and **EMT**. The goal is to visualize pattern formation, explore parameter effects, and analyze diffusion dynamics.

---

## Table of Contents
1. [Gray-Scott Model](#gray-scott-model)  
2. [Multi-fate Model](#multi-fate-model)  
3. [EMT Model](#emt-model)  
4. [Installation](#installation)  
5. [Usage](#usage)  
6. [Data & Results](#data--results)  
7. [Analysis & Insights](#analysis--insights)  
8. [Future Work](#future-work)

---

## Gray-Scott Model

**Description:**  
The Gray-Scott model simulates activator-inhibitor dynamics producing complex patterns.

### Model Dynamics
The system simulates the interaction of two chemical species, $u$ and $v$, where $u$ is replenished at a feed rate $F$ and $v$ is converted from $u$ in a self-replicating process:

$$
\begin{aligned}
\frac{\partial u}{\partial t} &= D_u \nabla^2 u - uv^2 + F(1-u) \\
\frac{\partial v}{\partial t} &= D_v \nabla^2 v + uv^2 - (F+k)v
\end{aligned}
$$

#### Parameter Definitions
* **$D_u, D_v$**: Diffusion coefficients for species $u$ and $v$.
* **$F$**: The "Feed" rate, determining how quickly $u$ is added to the system.
* **$k$**: The "Kill" rate, determining the rate at which $v$ is removed.

**Numerical Method:** Finite difference with explicit Euler integration.  
**Boundary Conditions:** Periodic   

**Sample Output:**  
![Gray–Scott Labyrinth Pattern](videos/Gray_Scott/Labirynths%20(F%20%3D%200.0545%2C%20k%3D0.062).gif)

See the folders Images and videos for more.

---

## Multi-fate Model

**Description:**  
The Multi-fate model simulates multiple cell fate decisions with diffusion-coupled reaction dynamics.

### Model Dynamics
The spatio-temporal evolution of species $A$ and $B$ is governed by the following reaction-diffusion system:

$$
\begin{aligned}
\frac{\partial A}{\partial t} &= D_u \nabla^2 A + f_A(A, B) \\
\frac{\partial B}{\partial t} &= D_v \nabla^2 B + f_B(A, B)
\end{aligned}
$$

#### Reaction Terms
The nonlinear reaction functions $f_A$ and $f_B$ describe the production and degradation rates:

$$
\begin{aligned}
f_A(A,B) &= \alpha + \frac{\beta t_A}{1 + t_A} - A \\
f_B(A,B) &= \alpha + \frac{\beta t_B}{1 + t_B} - B
\end{aligned}
$$

#### Auxiliary Variables
The Hill-like activation terms $t_A, t_B$ and the dimerization terms $A_2, B_2$ are defined as:

$$
t_A = (A_2)^n, \quad t_B = (B_2)^n
$$

$$
\begin{aligned}
A_2 &= \frac{2 A^2}{K + 4(A + B) + \sqrt{K^2 + 8 K_d (A + B)}} \\
B_2 &= \frac{2 B^2}{K + 4(A + B) + \sqrt{K^2 + 8 K_d (A + B)}}
\end{aligned}
$$

**Numerical Method:** Finite difference / Euler or other as implemented.  
**Boundary Conditions:** Periodic / Neumann

**Sample Output:**  
![MultiFate Radial Pattern](https://github.com/vajadiye-gif/-Reaction-Diffusion-Systems-as-Data-Generating-Dynamical-Systems/blob/main/videos/MultiFate/Radial%20Gradient.gif)

See the folders Images and videos for more.

---

## EMT Model

**Description:**  
The EMT (Epithelial-Mesenchymal Transition) model simulates diffusion-driven EMT patterns in cells.

### Model Dynamics
The spatio-temporal evolution of the regulatory species is governed by the following system of reaction-diffusion equations:

$$
\begin{aligned}
\frac{\partial [\text{mir200}]}{\partial t} &= G_{MIR200} \. h_{mir200}^{zeb} \. h_{mir200}^{sna} - [\text{mir200}]^1 ., \text{mir200}_0 - K_{MIR200} [\text{mir200}]^0 + D_{mir200} \nabla^2 [\text{mir200}] \\
\frac{\partial [m_{zeb}]}{\partial t} &= G_{MZEB} \. h_{mzeb}^{zeb} \. h_{mzeb}^{sna} - [\text{mir200}]^1 \. \text{mir200}_1 - K_{MZEB} [m_{zeb}] + D_{mzeb} \nabla^2 [m_{zeb}] \\
\frac{\partial [ZEB]}{\partial t} &= G_{ZEB} [\text{mir200}]^1 \. \text{mir200}_2 - K_{ZEB} [ZEB] + D_{ZEB} \nabla^2 [ZEB] \\
\frac{\partial [SNAIL]}{\partial t} &= G_{SNAIL} [\text{mir34}]^4 \. \text{mir34}_2 - K_{SNAIL} [SNAIL] + D_{SNAIL} \nabla^2 [SNAIL] \\
\frac{\partial [M_{SNAIL}]}{\partial t} &= G_{MSNAIL} \. h_{msnai} \, h_{msna}^{sna} - [\text{mir34}]^4 \. \text{mir34}_1 - K_{MSNAIL} [M_{SNAIL}] + D_{MSNAIL} \nabla^2 [M_{SNAIL}] \\
\frac{\partial [\text{mir34}]}{\partial t} &= G_{MIR34} \. h_{mir34}^{zeb} \. h_{mir34}^{sna} - [\text{mir34}]^4 \. \text{mir34}_0 - K_{MIR34} [\text{mir34}] + D_{mir34} \nabla^2 [\text{mir34}] \\
\frac{\partial [I]}{\partial t} &= 0
\end{aligned}
$$

#### Mathematical Definition of Phenotypes
We define the cellular phenotype $\mathcal{P}$ at any spatial coordinate $(i, j)$ based on the steady-state concentrations of the primary miR-200/ZEB loop:

$$
\mathcal{P}_{i,j} = 
\begin{cases} 
\text{Epithelial (E)} & \text{if } [\text{mir200}] \gg [\text{ZEB}] \\
\text{Hybrid (E/M)} & \text{if } [\text{mir200}] \approx [\text{ZEB}] \\
\text{Mesenchymal (M)} & \text{if } [\text{mir200}] \ll [\text{ZEB}] 
\end{cases}
$$

#### Phenotypic Signatures
The transition between states is characterized by the following molecular profiles:

$$
\begin{array}{|l|c|c|c|}
\hline
\textbf{Phenotype} & [\text{mir200}] & [\text{ZEB}] & \text{Biological State} \\
\hline
\color{blue}{\text{Epithelial (E)}} & \uparrow \text{ High} & \downarrow \text{ Low} & \text{Adherent / Stationary} \\
\hline
\color{green}{\text{Hybrid (E/M)}} & \sim \text{ Mid} & \sim \text{ Mid} & \text{Collective Migration} \\
\hline
\color{red}{\text{Mesenchymal (M)}} & \downarrow \text{ Low} & \uparrow \text{ High} & \text{Invasive / Solitary} \\
\hline
\end{array}
$$

**Numerical Method:** Finite difference / Euler.  
**Boundary Conditions:** Neumann

**Sample Output:**  
![EMT SNAIL Pattern](https://github.com/vajadiye-gif/-Reaction-Diffusion-Systems-as-Data-Generating-Dynamical-Systems/blob/main/videos/EMT/SNAIL.gif)

See the folders Images and videos for more.

---

## Installation

pip install -r requirements.txt

## Usage

Open the corresponding notebook in notebooks/ (in Google Colab for interactivity)

Select the species using ToggleButtons (e.g., u, v, SNAIL, ZEB ,etc.)
Run the notebook → visualize patterns
Save results as arrays, images, or gifs if needed.

## Data & Results

Output videos: videos/

Sample images: Images/

Interactive widgets allow parameter exploration in real-time.

**Note on Notebooks:**
Due to interactive widgets and long-running simulations, the notebooks are also provided as exported PDFs for static viewing on GitHub.
The original interactive versions are intended to be run in Google Colab.

## Analysis & Insights

Gray-Scott: Increasing F accelerates pattern formation.
Multi-fate: Different diffusion constants produce distinct fate domains.
EMT: Varying reaction rates changes cluster morphology.

## Future Work

Add more reaction-diffusion models.
Explore 3D pattern formation.
Optimize simulation speed for larger grids.
Compare parameter sensitivity across models.
