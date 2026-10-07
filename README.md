# Dynamics, bifurcations and extensions of models of unified growth

## 🚀 Description 

Unified growth theory (UGT) is family of models proposed in [[Galor 2000]](https://www.aeaweb.org/articles?id=10.1257/aer.90.4.806) and used by macroeconomists to explain the three main stages of population growth and technological progress throughout civilisations.
At its core, UGT derives the laws of motions by maximizing the utility of each agent in the adult population.
Each agent choses the optimal allocation of resources between consumption and child-rearing activities, depending on their income level and the state of the economy.

Such problem yields, for each generation t>0, an optimal number of children an agent will have and the level of education they bestow upon them, raising the human capital in the economy at t+1.
From Galor's analysis, this gives rise to a four-dimensional, piecewise continuous and non-smooth iterated map.
Later efforts [[Lagerlöf 2006]](https://www.aeaweb.org/articles?id=10.1257/aer.90.4.806) have tested this model using numerical simulations, however they were restricted to a single set of parameters and a single initial state for a starting economy.

### ✏️  Outline

__Collaborators:__ Samuel Bolduc St-Aubin (U. of Auckland, New Zealand), Sam Doak (U. of Auckland, New Zealand) and Greta Meggiorini (U. of Auckland, New Zealand). 

We certify Lagerlöf's numerical findings and expand upon them by putting UGT models through the rigorous machinery of dynamical systems theory and bifurcation analysis.
Among the numerous findings, we hereby list the most important ones:

- the family of models proposed by Galor is actually three-dimensional, not four dimensional (education is a choice variable rather than an endogenous one);
- the [Malthusian trap](https://en.wikipedia.org/wiki/Malthusianism#Theory_of_breakout_via_technology) is not a fixed point for the family of models proposed by Galor, making claims such as the following one incorrect
> (From _"UNIFIED GROWTH THEORY: ROOTS OF GROWTH AND INEQUALITY IN THE WEALTH OF NATIONS", page 9, third paragraph_) ``Unified Growth Theory uncovers the central forces that gave rise to the Malthusian trap, unraveling the mechanisms that ultimately enabled humanity to escape its gravitational pull''
- Lagerlöf's modeling choices yield a system with a substantial amount of economies that undergo population collapse within one or two generations, despite [[Lagerlöf 2006]](https://www.aeaweb.org/articles?id=10.1257/aer.90.4.806) showing a single simulation for which the economy reaches the sustained growth regime;
- Lagerlöf's specified parameter values are not robust, as the system sits at a border-collision bifurcation (see [[Simpson 2022](https://www.tandfonline.com/doi/full/10.1080/10236198.2023.2265495)]);
- the steady-state in the modern growth regime is explained by looking at the large technology asymptotic limit of Lagerlöf's system; in this regime the dynamics reduces to a two-dimensional state space whose partition into a 6-subsets covering, and transitions therein, induces a directed graph that explains the different endstates of any economy that survives population collapse.

### 💡 What does each experiment do?

#### fig_01.m
Computes the illustrative trajectory of the GWL model replicating the calibration done in [[Lagerlöf 2006]](https://www.aeaweb.org/articles?id=10.1257/aer.90.4.806). 

#### fig_02.m
Computes an illustrative trajectory of the GWL model with a slightly perturbed value of $L_0$, showing the population collapse following the event $L_1 = 0$.

#### fig_03.m
Classifies the initial conditions in the two-dimensional $(L_0, A_0)$-slice, at chosen $G_0$, of the full state space $\mathbb{R}^3_{\geq 0}$ of the GWL model.

#### fig_04.m
Plots the boundaries of the switching manifolds of the reduced state space at selected $G_0$ slices obtained via `build_M1_M2.m`. 

#### fig_07.m
Computes three illustrative trajectories of economies in the sustained growth regime at varying values of the preference parameter $\gamma$.

#### fig_08a.m
Computes the bifurcation diagram of the sustained growth steady-state of the asymptotic reduction in the preference parameter $\gamma$. 

#### fig_08b.m
Computes and plots the three trajectories of `fig_07.m` but in the partitioned state space in the asymptotic sustained growth regime.

#### fig_09.m
Computes an illustrative trajectories of the sGW model.

#### fig_10.m
Similarly to `fig_03.m`, it classifies the initial conditions in the two-dimensional $(L_0,A_0)$-slice, at chosen $G_0$, of the full state space $\mathbb{R}^3_{\geq 0}$ of the sGW model.

#### build_M1_M2.jl
Julia script that builds and exports the boundaries $M_1$ and $M_2$ of the regions $U_1$ and $U_2$ of population collapse (blue and pink region). The results are found in the [boundaries](./boundaries/) subfolder.
