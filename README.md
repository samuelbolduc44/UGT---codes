# Dynamics, bifurcations and extensions of models of unified growth

## 🚀 Description 

Unified growth theory (UGT) is family of models proposed in [[Galor 2000]](https://www.aeaweb.org/articles?id=10.1257/aer.90.4.806) and used by macroeconomists to explain the three main stages of population growth and technological progress throughout civilisations.
At its core, UGT derives the equations of motions by optimising the utility of each agent in the working population.
Each agent choses the optimal allocation of resources between consumption and child-rearing activities, depending on their income level and the endogenous state of the economy.

Such problem yields, for each generation t>0, an optimal number of children an agent will have and the level of education they bestow upon them, raising the human capital in the economy at t+1.
From Galor's analysis, this gives rise to a four-dimensional, piece-wise continuous (i.e. non-smooth) iterated map.
Later efforts [[Lagerlöf 2006]](https://www.aeaweb.org/articles?id=10.1257/aer.90.4.806) have tested this model using numerical simulations, however they stopped at one set of parameters and a single starting economy.

### ✏️  Outline

__Collaborators:__ Samuel Bolduc St-Aubin (U. of Auckland, New Zealand), Sam Doak (U. of Auckland, New Zealand) and Greta Meggiorini (U. of Auckland, New Zealand). 

We certify Lagerlöf's numerical findings and expand upon them by putting UGT models through the rigorous machinery of dynamical systems theory and bifurcation analysis.
Among the numerous findings, we hereby list the most important ones:

- the family of models proposed by Galor is actually three-dimensional, not four dimensional (education is a choice variable rather than an endogenous one);
- the [Malthusian trap](https://en.wikipedia.org/wiki/Malthusianism#Theory_of_breakout_via_technology) is not a fixed point for the family of models proposed by Galor, making claims such as the following one incorrect
> (From _"UNIFIED GROWTH THEORY: ROOTS OF GROWTH AND INEQUALITY IN THE WEALTH OF NATIONS", page 9, third paragraph_) ``Unified Growth Theory uncovers the central forces that gave rise to the Malthusian trap, unraveling the mechanisms that ultimately enabled humanity to escape its gravitational pull''
- Lagerlöf's modelling choices yield a system with a substantial amount of economies that undergo population collapse within one or two generations, despite [[Lagerlöf 2006]](https://www.aeaweb.org/articles?id=10.1257/aer.90.4.806) showing a single simulation for which the economy reaching post-modern growth;
- Lagerlöf's specified parameter values are not robust, as the system sits at a border-collision bifurcation (see [[Simpson 2022](https://www.tandfonline.com/doi/full/10.1080/10236198.2023.2265495)]);
- the post-modern growth regime is explained by looking at the large technology asymptotic limit of Lagerlöf's system; in this regime the dynamics reduces to a two-dimensional state space whose partition into a 6-subsets covering, and transitions therein, induces a directed graph that explains the different endstates of any economy that survives population collapse.

### 💡 What does each experiment do?

#### FIG1.m
An illustrative trajectory of the GWL model replicating the calibration done by Lagerlöf (2008). The economy, represented by the endogenous state vector $X_t = (G_t, A_t, L_t)$ (panels (a)-(c)) and the endogenous variables $e_t$, $z_t$, $h_t$ (panels (d)-(f)), undergoes the three stages predicted by unified growth theory. 

#### FIG2.m
An illustrative trajectory of the GWL model with a slightly perturbed value of $L_0$, showing the population collapse following the event $L_1 = 0$.

#### varyinggamma.m
Solve asymptotic system. Used to make FIG7

#### phaseportraitandtraj.m
Panel (b) of FIG8. 

#### bifdiagram_paper.m
Panel (a) of FIG8.

#### Counting_itr_before_collapse.m
Helper function for producing FIG 3 and FIG 4. Simply adjust the range of $A_0$, $L_0$, and the value for $G_0$.

Specifically, it computes the classification of initial conditions in the two-dimensional $(L_0, A_0)$-slice, at chosen $G_0$, of the full state space $\mathbb{R}^3_{\geq 0}$.

- **White region:** Initial conditions that do not undergo the collapse set over the simulated horizon.
- **Blue region:** Initial conditions that undergo population collapse after one generation. 
- **Pink region:** Initial conditions that undergo population collapse after two generations.

The code uses a parallel loop for efficiency.

#### makeboundaries.m
Matlab script that plots the boundaries obtained via build_M1_M2.m resulting in Figure 4. 

#### build_M1_M2.m
Julia script that builds, plots and exports (in .mat format) the boundaries $\partial U_1$ and $\partial U_2$ of the regions $U_1$ and $U_2$ of population collapse (blue and pink region) computed empirically by 'Counting_itr_before_collapse.m'
