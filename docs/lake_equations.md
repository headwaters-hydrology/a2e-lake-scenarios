# Lake equations for changes in input concentrations
These are the simplified lake equations solved for the change ratio as defined as the new input concentrations (from reductions to mitigations) to the original input concentrations. These are derived from Abell et al 2019 and Abell et al 2020.

### Definition of the change ratio

$$
r_{in()} = \frac{C_{in()}^{scenario}}{C_{in()}^{current}}
$$

where...

$C_{in()}^{scenario}$ is the inflow concentration of either TN or TP in the scenario (with the mitigations) and $C_{in()}^{current}$ is the inflow concentration of either TN or TP in the current conditions.

### Application of the lake change ratios
Once a calculation of $r_{lake()}$ has been made, then you need to multiply it by the current lake concentration ($C_{lake()}^{current}$) to get the scenario lake concentration ($C_{lake()}^{scenario}$):

$$
C_{lake()}^{scenario} = r_{lake()} C_{lake()}^{current}
$$

### Examples
Each of the following equations will have examples using actual data to compare the results of the simplified equations to the original Abell equations.

The source data variables include:

$$
C_{in(P)}^{current} = 100 \space \frac{mg}{m^3}
$$

$$
C_{in(P)}^{scenario} = 40 \space \frac{mg}{m^3}
$$

$$
C_{in(N)}^{current} = 4000 \space \frac{mg}{m^3}
$$

$$
C_{in(N)}^{scenario} = 800 \space \frac{mg}{m^3}
$$

$$
z_{max} = 20 \space m
$$

$$
\tau = 10 \space years
$$

subsequently...

$$
r_{in(P)} = \frac{0.04}{0.1} = 0.4
$$

$$
r_{in(N)} = \frac{0.8}{4} = 0.2
$$

## Change ratios

### TP
#### Original equation
When $\hspace{2mm} z_{max} \gt 7.5$:

$$
\log_{10} C_{lake(P)} = \frac{\log_{10} C_{in(P)}}{1 + 0.44\tau^{0.13}}
$$

else:

$$
\log_{10} C_{lake(P)} = \log_{10} C_{in(P)}
$$

#### Simplified equations
When $\hspace{2mm} z_{max} \gt 7.5$:

$$
r_{lake(P)} = r_{in(P)}^{1/b}
$$

where...

$$
b = 1 + 0.44 \tau^{0.13}
$$

$r_{in(P)}$ is the TP concentration change ratio in the inflows to the lake, $r_{lake(P)}$ is the in-lake TP concentration change ratio, $z_{max}$ is the max depth of the lake in meters, and $\tau$ is the residence time of the lake in years.

Else when $\hspace{2mm} z_{max} \lt= 7.5$:

$$
r_{lake(P)} = r_{in(P)}
$$

#### Example

$$
b = 1 + 0.44 \times 10^{0.13} = 1.5935
$$

$$
\log_{10} C_{lake(P)}^{current} = \frac{\log_{10} 100}{1.5935} = 1.255
$$

$$
C_{lake(P)}^{current} = 10^{1.255} = 18
$$

$$
\log_{10} C_{lake(P)}^{scenario} = \frac{\log_{10} 40}{1.5935} = 1.005
$$

$$
C_{lake(P)}^{scenario} = 10^{1.005} = 10.124
$$

$$
\frac{C_{lake(P)}^{scenario}}{C_{lake(P)}^{current}}= \frac{10.124}{18} = 0.5625
$$

$$
r_{lake(P)} = 0.4^{1/1.5935} = 0.562
$$

### TN
#### Original equation

$$
\log_{10} C_{lake(P)} = 1.6 + 0.54 \log_{10} C_{in(N)} - 0.41 \log_{10} z_{max}
$$

#### Simplified equation

$$
r_{lake(N)} = r_{in(N)}^{0.54}
$$

#### Example

$$
\log_{10} C_{lake(N)}^{current} = 1.6 + 0.54 \log_{10} 4000 - 0.41 \log_{10} 20  = 3.012
$$

$$
C_{lake(N)}^{current} = 10^{3.012} = 1028
$$

$$
\log_{10} C_{lake(N)}^{scenario} = 1.6 + 0.54 \log_{10} 800 - 0.41 \log_{10} 20  = 2.634
$$

$$
C_{lake(N)}^{scenario} = 10^{2.634} = 430.5
$$

$$
\frac{C_{lake(N)}^{scenario}}{C_{lake(N)}^{current}} = \frac{430.5}{1028} = 0.419
$$

$$
r_{lake(N)} = 0.2^{0.54} = 0.419
$$

### Chla
#### Original equation

$$
\log_{10} C_{lake(Chla)} = -1.8 + 0.7 \log_{10} C_{lake(N)} + 0.55 \log_{10} C_{lake(P)}
$$

#### Simplified equation

$$
r_{lake(Chla)} = r_{lake(N)}^{0.7} r_{lake(P)}^{0.55}
$$

#### Example

$$
\log_{10} C_{lake(Chla)}^{current} = -1.8 + 0.7 \times 3.012 + 0.55 \times 1.255 = 0.999
$$

$$
C_{lake(Chla)}^{current} = 10^{0.999} = 9.977
$$

$$
\log_{10} C_{lake(Chla)}^{scenario} = -1.8 + 0.7 \times 2.634 + 0.55 \times 1.005 = 0.5965
$$

$$
C_{lake(Chla)}^{scenario} = 10^{0.5965} = 3.949
$$

$$
\frac{C_{lake(Chla)}^{scenario}}{C_{lake(Chla)}^{current}} = \frac{3.949}{9.977} = 0.396
$$

$$
r_{lake(Chla)} = 0.419^{0.7} 0.562^{0.55} = 0.396
$$

### Secchi
#### Original equation
When $\hspace{2mm} z_{max} \geq 20$:

$$
D_{lake(Secchi)}^{0.5} = 3.46 - 1.53 \log_{10} C_{lake(Chla)}
$$

else:

$$
D_{lake(Secchi)}^{0.5} = 3.46 - 0.74 \log_{10} C_{lake(Chla)} - 0.35 \log_{10} \frac{Fetch U^2}{z_{max}}
$$

#### Simplified equations

When $\hspace{2mm} z_{max} \geq 20$:

$$
D_{lake(Secchi)}^{scenario} = (\log_{10} (r_{lake(Chla)}^{-1.53}) + (D_{lake(Secchi)}^{current})^{0.5})^2
$$

Else when $\hspace{2mm} z_{max} \lt 20$:

$$
D_{lake(Secchi)}^{scenario} = (\log_{10} (r_{lake(Chla)}^{-0.74}) + (D_{lake(Secchi)}^{current})^{0.5})^2
$$

#### Example

$$
\sqrt D_{lake(Secchi)}^{current} = 3.46 - 1.53 \times 0.999 = 1.931
$$

$$
D_{lake(Secchi)}^{current} = 1.931^2 = 3.729
$$

$$
\sqrt D_{lake(Secchi)}^{scenario} = 3.46 - 1.53 \times 0.5965 = 2.547
$$

$$
D_{lake(Secchi)}^{scenario} = 2.547^2 = 6.487
$$

$$
\frac{D_{lake(Secchi)}^{scenario}}{D_{lake(Secchi)}^{current}} = \frac{6.487}{3.729} = 1.734
$$

$$
D_{lake(Secchi)}^{scenario} = (\log_{10} (0.396^{-1.53}) + 3.729^{0.5})^2 = 6.487
$$


## Waikato and Canterbury are special
### TP Waikato

$$
b = 1 + \tau^{0.5}
$$

#### Original equation

$$
\log_{10} C_{lake(P)} = 0.9217 + 0.6172 \frac{\log_{10} C_{in(P)}}{b}
$$

#### Simplified equation

$$
r_{lake(P)} = r_{in(P)}^{0.6172/b}
$$

### TP Canterbury

When $\hspace{2mm} z_{max} \gt 7.5$:

$$
b = 1 + 0.91888 t^{0.0205}
$$

else:

$$
b = 1 + 0.09288 t^{0.0205}
$$

#### Original equation

$$
\log_{10} C_{lake(P)} = \frac{\log_{10} C_{in(P)}}{b}
$$

#### Simplified equation

$$
r_{lake(P)} = r_{in(P)}^{1/b}
$$

### TN Waikato

$$
b = 1 + \tau^{0.5}
$$

#### Original equation

$$
\log_{10} C_{lake(N)} = 2.3969 + 0.3564 \frac{\log_{10} C_{in(N)}}{b}
$$

#### Simplified equation

$$
r_{lake(N)} = r_{in(N)}^{0.3564/b}
$$



























