# Codes Description

## FIG1.m
An illustrative trajectory of the GWL model replicating the calibration done by Lagerlöf (2008). The economy, represented by the endogenous state vector $X_t = (G_t, A_t, L_t)$ (panels (a)-(c)) and the endogenous variables $e_t$, $z_t$, $h_t$ (panels (d)-(f)), undergoes the three stages predicted by unified growth theory. 

## FIG2.m
An illustrative trajectory of the GWL model with a slightly perturbed value of $L_0$, showing the population collapse following the event $L_1 = 0$.

## Counting_itr_before_collapse.m
Helper function for producing FIG 3 and FIG 4. Simply adjust the range of $A_0$, $L_0$, and the value for $G_0$.

Specifically, it computes the classification of initial conditions in the two-dimensional $(L_0, A_0)$-slice, at chosen $G_0$, of the full state space $\mathbb{R}^3_{\geq 0}$.

- **White region:** Initial conditions that do not undergo the collapse set over the simulated horizon.
- **Blue region:** Initial conditions that undergo population collapse after one generation. 
- **Pink region:** Initial conditions that undergo population collapse after two generations.

The code uses a parallel loop for efficiency.
