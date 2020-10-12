# Empires

A simulation of human interactions at a large scale, using python numpy and matplotlib.
It is a simple reaction–diffusion system based on darwinian dynamics.

![](documentation/example.gif)


# How does-it work ?

This script solves the system:

<a href="https://www.codecogs.com/eqnedit.php?latex=\left\{&space;\begin{array}{ll}&space;\partial_tN_i&space;=&space;G_i(N,&space;\rho,&space;v,&space;u)|_{v=u_i}&space;N_i&space;&plus;&space;D_i\nabla&space;e^{-G_i(N,&space;\rho,u)}\nabla&space;N_i&space;&plus;&space;d_i&space;N_i&space;\Delta&space;(\sum_j&space;\rho_j)&space;\\&space;\partial_t\rho_j&space;=&space;r&space;\rho_j&space;(1-\frac{\rho_j}{K_j})&space;-&space;\sum_i{&space;\Phi_i(\rho_j)&space;N_i&space;}&space;\\&space;\partial_t&space;u_i&space;=&space;\partial_v&space;G_i(N,&space;\rho,&space;v,&space;u)|_{v=u_i}&space;\end{array}&space;\right." target="_blank"><img src="https://latex.codecogs.com/gif.latex?\left\{&space;\begin{array}{ll}&space;\partial_tN_i&space;=&space;G_i(N,&space;\rho,&space;v,&space;u)|_{v=u_i}&space;N_i&space;&plus;&space;D_i\nabla&space;e^{-G_i(N,&space;\rho,u)}\nabla&space;N_i&space;&plus;&space;d_i&space;N_i&space;\Delta&space;(\sum_j&space;\rho_j)&space;\\&space;\partial_t\rho_j&space;=&space;r&space;\rho_j&space;(1-\frac{\rho_j}{K_j})&space;-&space;\sum_i{&space;\Phi_i(\rho_j)&space;N_i&space;}&space;\\&space;\partial_t&space;u_i&space;=&space;\partial_v&space;G_i(N,&space;\rho,&space;v,&space;u)|_{v=u_i}&space;\end{array}&space;\right." title="\left\{ \begin{array}{ll} \partial_tN_i = G_i(N, \rho, v, u)|_{v=u_i} N_i + D_i\nabla e^{-G_i(N, \rho,u)}\nabla N_i + d_i N_i \Delta (\sum_j \rho_j) \\ \partial_t\rho_j = r \rho_j (1-\frac{\rho_j}{K_j}) - \sum_i{ \Phi_i(\rho_j) N_i } \\ \partial_t u_i = \partial_v G_i(N, \rho, v, u)|_{v=u_i} \end{array} \right." /></a>

## First equation: population evolution

The first equation of the system represents the evolution of each population:
+ _G_ is the reproduction rate of the population in a space cell (also called numerical response, or G-function).
+ _D_ is the natural diffusion coefficient, the exponential term is here to model the fact that emigration is more important from cells where the reproduction rate is negative
+ _d_ is the drift coefficient, which models the capacity of planned expension (the higher it is, the more the pop will be attracted to places where resources are numerous)

### The G-function
The G-function has the form <a href="https://www.codecogs.com/eqnedit.php?latex=G_i(N,&space;\rho,&space;v,&space;u)&space;=&space;G^1_i(\rho)&space;&plus;&space;G^2_i(N,v,u)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?G_i(N,&space;\rho,&space;v,&space;u)&space;=&space;G^1_i(\rho)&space;&plus;&space;G^2_i(N,v,u)" title="G_i(N, \rho, v, u) = G^1_i(\rho) + G^2_i(N,v,u)" /></a>, with the first term modelling reproduction without any external perturbations, the second the influence of other populations.

For the first term, two main forms are used in population dynamics:
+ Laissez-faire: <a href="https://www.codecogs.com/eqnedit.php?latex=G^1_i(\rho)&space;=&space;\frac{b}{c}&space;\sum_j&space;{\Phi_i(\rho_j)}&space;-&space;\gamma" target="_blank"><img src="https://latex.codecogs.com/gif.latex?G^1_i(\rho)&space;=&space;\frac{b}{c}&space;\sum_j&space;{\Phi_i(\rho_j)}&space;-&space;\gamma" title="G^1_i(\rho) = \frac{b}{c} \sum_j {\Phi_i(\rho_j)} - \gamma" /></a>

+ Leslie: <a href="https://www.codecogs.com/eqnedit.php?latex=G^1_i(\rho)&space;=&space;s(1-\frac{N_i}{h&space;\sum_j{\rho_j}})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?G^1_i(\rho)&space;=&space;s(1-\frac{N_i}{h&space;\sum_j{\rho_j}})" title="G^1_i(\rho) = s(1-\frac{N_i}{h \sum_j{\rho_j}})" /></a>

The Laissez-faire model is used here.

The second term depends from a vector _u_ which represents the war strategy of each population which evolves according to darwinian dynamics (third equation of the main system). And we have (see Vincent, Thomas & all):

+ <a href="https://www.codecogs.com/eqnedit.php?latex=G^2_i(N,&space;v,&space;u)&space;=&space;-&space;\frac{1}{K_N}&space;\sum_k&space;e^{-\frac{(v-u_k)^2}{2&space;\sigma_\alpha}}&space;N_k" target="_blank"><img src="https://latex.codecogs.com/gif.latex?G^2_i(N,&space;v,&space;u)&space;=&space;-&space;\frac{1}{K_N}&space;\sum_k&space;e^{-\frac{(v-u_k)^2}{2&space;\sigma_\alpha}}&space;N_k" title="G^2_i(N, v, u) = - \frac{1}{K_N} \sum_k e^{-\frac{(v-u_k)^2}{2 \sigma_\alpha}} N_k" /></a>




## Second equation: resources evolution

The resources have a natural growth rate r, are limited by a carrying capacity K and are consumed by pops.

### The functional response

Consumption is represented by functional response. In population dynamics, thoses functional responses were described by Holling and can have three forms:

+ Holling type I: <a href="https://www.codecogs.com/eqnedit.php?latex=\Phi_I(\rho)&space;=&space;\left\{&space;\begin{array}{ll}&space;\frac{c}{2\rho_0}&space;&&space;\mbox{if&space;}&space;\rho&space;<&space;2&space;\rho_0&space;\\&space;c&space;&&space;\mbox{else}&space;\end{array}&space;\right." target="_blank"><img src="https://latex.codecogs.com/gif.latex?\Phi_I(\rho)&space;=&space;\left\{&space;\begin{array}{ll}&space;\frac{c}{2\rho_0}&space;&&space;\mbox{if&space;}&space;\rho&space;<&space;2&space;\rho_0&space;\\&space;c&space;&&space;\mbox{else}&space;\end{array}&space;\right." title="\Phi_I(\rho) = \left\{ \begin{array}{ll} \frac{c\rho}{2\rho_0} & \mbox{if } \rho < 2 \rho_0 \\ c & \mbox{else} \end{array} \right." /></a>

+ Holling type II: <a href="https://www.codecogs.com/eqnedit.php?latex=\Phi_{II}(\rho)&space;=&space;\frac{c\rho}{\rho_0&space;&plus;&space;\rho}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\Phi_{II}(\rho)&space;=&space;\frac{c\rho}{\rho_0&space;&plus;&space;\rho}" title="\Phi_{II}(\rho) = \frac{c\rho}{\rho_0 + \rho}" /></a>

+ Holling type III: <a href="https://www.codecogs.com/eqnedit.php?latex=\Phi_{II}(\rho)&space;=&space;\frac{c\rho^2}{\rho_0^2&space;&plus;&space;\rho^2}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\Phi_{III}(\rho)&space;=&space;\frac{c\rho^2}{\rho_0^2&space;&plus;&space;\rho^2}" title="\Phi_{II}(\rho) = \frac{c\rho^2}{\rho_0^2 + \rho^2}" /></a>

For the moment, only Holling III is implemented by default.

## Third equation: darwinian dynamics

The third equation represents darwinian evolution of war strategies. It will probably be replaced by a more detailed game theory model in the future.



# To-do
- [ ] Use cuda and cupy for better perf
- [ ] Try different solving methods to gain perfs
- [x] Tune parameters
- [x] Use more coherent game theory models
- [ ] Generate resources from maps
- [ ] Make a clear documentation on the maths behind the model

# Sources

On taxis and kinesis equations:
+ Arthur, R. M. « SPECIES PACKING, AND WHAT COMPETITION MINIMIZES ». _Proceedings of the National Academy of Sciences_, vol. 64, no 4, décembre 1969, p. 1369‑71. DOI.org (Crossref), doi:10.1073/pnas.64.4.1369.
+ Averill, Isabel, et al. « On Several Conjectures from Evolution of Dispersal ». _Journal of Biological Dynamics_, vol. 6, no 2, mars 2012, p. 117‑30. DOI.org (Crossref), doi:10.1080/17513758.2010.529169.
+ Gorban, A. N., et N. Çabukoǧlu. « Basic Model of Purposeful Kinesis ». _Ecological Complexity_, vol. 33, janvier 2018, p. 75‑83. DOI.org (Crossref), doi:10.1016/j.ecocom.2018.01.002.
+ Hillen, T., et K. J. Painter. « A User’s Guide to PDE Models for Chemotaxis ». _Journal of Mathematical Biology_, vol. 58, no 1‑2, janvier 2009, p. 183‑217. DOI.org (Crossref), doi:10.1007/s00285-008-0201-3.
+ Leeuwen, E. & Jansen, Vincent & Bright, Paul. (2007). How population dynamics shape the functional response in a one-predator-two-prey system. Ecology. 88. 1571-81. 10.1890/06-1335.
+ Vincent, Thomas. (2004). The G-function method for analyzing Darwinian dynamics. International Game Theory Review. 6. 10.1142/S0219198904000083.

On evolutionary game theory:
+ Tembine, Hamidou, et al. « Evolutionary Games in Wireless Networks ». _IEEE Transactions on Systems, Man, and Cybernetics, Part B (Cybernetics)_, vol. 40, no 3, juin 2010, p. 634‑46. DOI.org (Crossref), doi:10.1109/TSMCB.2009.2034631.
+ Cressman R., Apaloo J. (2016) Evolutionary Game Theory. In: Basar T., Zaccour G. (eds) Handbook of Dynamic Game Theory. Springer, Cham. https://doi.org/10.1007/978-3-319-27335-8_6-1
+ Perc, Matjaž, et al. « Statistical Physics of Human Cooperation ». _Physics Reports_, vol. 687, mai 2017, p. 1‑51. DOI.org (Crossref), doi:10.1016/j.physrep.2017.05.004.
+ Sandholm, William H. « Evolutionary Game Theory ». _Encyclopedia of Complexity and Systems Science_, édité par Robert A. Meyers, Springer Berlin Heidelberg, 2017, p. 1‑38. DOI.org (Crossref), doi:10.1007/978-3-642-27737-5_188-3.
