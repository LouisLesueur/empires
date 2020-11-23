---
title: Empires documentation
---

# Introduction: real modeling or just making equations ?

Clyodynamics is a research field essentially carried by Peter Turchin. Its goal is to bring mathematical modelisation into history research. The most complete and detailed presentation of Turchin theory is probably in his book *Historical Dynamics: Why States Rise and Fall*[^1]. Turchin seems to find interesting results, and especially his 'secular cycles'. Today, cliodynamics is still a new and growing field of research. But is Turchin analysis relevent ?

As recently underlined by A.Maini [^3], Turchin models suffer from a lack of strong theoretical background. Especially all the logistic growth models used in its population dynamics equations relies on many parameters that are just impossible to estimate. And this is mostly due to the fact that population dynamics is not a fully rigorous theory, as many of the equations used in biology are built to fit with the observation and not from general principles.

Despite these observations, cliodynamics stays a very interesting field, and equations modelling population evolution are not its only elements. But it is important to stay rigorous and to not get carried away by results that just seems to be right.

For all those reasons, the model presented in this paper, even if a part of it relies on scientific papers, need to be seen only as a toy model without any pretention.

![Example of simulation](example.gif)


# Space partition

## Basic entities

+ The full domain (map) of the simulation is noted $\Omega$. It is divided between water $\mathcal{W}$ and land $\mathcal{L}$.
+ The resource function $R(x,y,t)$ is defined on $\mathcal{L}$ and represents the spatial repartition of resources on the map over time.
+ A city is a point $c \in \mathcal{L}$. It has a zone of influence $I(c)$, centered around it, in which it can collect resources. For a given map, the set of cities is notes $\mathcal{C}(\Omega)$
+ The population function $N(c)$ is defined on $\mathcal{C}$ and represents the population of each city.

## Politic entities

+ A city policy is axed on three possible actions:
    + build a 'road' and a colony $c'$ in its expansion area $E(c)$.
    + disappear (because of war or collapse, see next sections)
    + build a 'road' to another city
+ A state is a connected graph of cities, linked by roads. It evolves over time, starting by one unique city (called a city-state) which can found colonies or conquer other states cities to expand. A state policy is defined by:
    + Collecting taxes
    + Spending public money in various domains
    + Defining its attitude (ie the attitude of its cities) towards other states

# Demography

The evolution of population, public resource ($R_p$) and resources of a city $c$ is driven by a generalized prey-predator equation of the type:

$$
\left\{
    \begin{array}{ll}
        \frac{dN}{dt} = \chi f(N,R) - \delta(N)N & \text{on } c\\
        \frac{dR_p}{dt} = \alpha(r(R)R - f(N,R)) & \text{on } c \\
        \frac{dR}{dt} = (1-\alpha)(r(R)R - f(N,R)) & \text{on } I(c)
    \end{array}
\right.
$$

On the rest of land, resources evolve according to $\frac{dR}{dt} = r(R)R$

where:

+ $r$ is the resource growth in the absence of predators
+ $\delta$ is the population decline in the absence of resource
+ $f$ is the population functional response (consumption function)
+ $\chi$ is the conversion rate of resources into people
+ $\alpha$ is the taxation rate

Several population dynamics models can be found in Turchin book on population dynmaics[^2]. For the moment only the one called Bazykin model by Turchin is implemented, but others will be added in the future.

All the parameters of the system are defined on the domain (ie depends on the position of the city), in function of a 'resource map' and the topograhic map which are needed in the code input.

# Expansion

In the first versions of the code, migration was simulated by directly adding a migration term in the eqation. Even if this kind of modelization is well documented[^4], it was abandonned. Indeed, classical migration terms are taken from diffusion equations. If for animal migration it can have an interest, it is not adapted to model human behaviour, and in particular the capacity of humans to group in a city which can exploit a large amount of land to feed its inhabitants. Adopting the city model is also better to properly define states and policies.

So, for a given city $c$, at each time step, a random point in $e(c)$ is sampled. Then, the probabilty to establish a colony on this point is computed, and if it is superior to a given threshold (which is decided by the state as a policy), the new colony is founded. For the moment, evry probabilistic step is uniform and a city is always founded. But this will change in the future.

# Politics

### State evolution

At $t=0$, the map is composed of city-states. At each time step, they can expand by creating colonies according to the rules mentionned above. An initial city-states and all its descendents form the state.

### Migration and trade

Inside a state, migrating and sending resources from a city to another is possible. Proper rules to compute the probability to do this need to be established.

### Rebellion
The satisfaction of population is defined on cities as:
$$
S(c) = R(I(c)) - R_m
$$
where $R_m$ is the minimal resources needed by the city population to 'be happy'. When $S(c) <0$, $c$ goes into rebellion and becomes a new city-state which forms the basis of a new state.

### War and peace
When a state discovers the existence of another states (which happens when cities are trying to expand onto the other state's land), it can do several actions:
+ Attack (in this case, the state with higher military budget has a higher probabilty to win)
+ Become ally and build a road to allow migrations and trade.



[^1]: Turchin, Peter. Historical Dynamics: Why States Rise and Fall: Why States Rise and Fall. Princeton University Press, 2003. JSTOR, www.jstor.org/stable/j.ctt1wf4d45. Accessed 23 Nov. 2020.
[^2]: Turchin, P.B.. (2003). Complex population dynamics: A theoretical/empirical synthesis. 10.1515/9781400847280.
[^3]: Maini, A. On Historical Dynamics by P. Turchin. Biophys Econ Sust 5, 3 (2020). https://doi.org/10.1007/s41247-019-0063-x
[^4]: Cui, Shangbin & Bai, Meng. (2014). Mathematical analysis of population migration and its effects to spread of epidemics. Discrete and Continuous Dynamical Systems - Series B. 20. 10.3934/dcdsb.2015.20.2819.
