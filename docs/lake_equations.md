## Lake equations for changes in input concentrations
These are the simplified lake equations solved for the change ratio as defined as the new input concentrations (from reductions to to mitigations) to the original input concentrations. These are derived from Abell et al 2019 and Abell et al 2020.

### TP
#### Current conditions
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

#### "Reference" conditions
When $\hspace{2mm} z_{max} \gt 7.5$:

$$
r_{lake(P)} = r_{in(P)}^{1/b_{ref}}
$$

where...

$$
b_{ref} = 1 + 0.27 \tau^{0.29}
$$

Else when $\hspace{2mm} z_{max} \lt= 7.5$:

$$
r_{lake(P)} = r_{in(P)}
$$

### TN
#### Current conditions

$$
r_{lake(N)} = r_{in(N)}^{0.54}
$$

#### "Reference" conditions

$$
r_{lake(N)} = r_{in(N)}^{0.81}
$$

### Chla
#### Current conditions

$$
r_{lake(Chla)} = r_{lake(N)}^{0.7} r_{lake(P)}^{0.55}
$$

#### "Reference" conditions

$$
r_{lake(Chla)} = r_{lake(N)}^{0.65} r_{lake(P)}^{0.59}
$$

### Secchi
#### Current conditions
When $\hspace{2mm} z_{max} \gt 20$:

$$
r_{lake(Secchi)} = r_{lake(Chla)}^{0.9}
$$

Else when $\hspace{2mm} z_{max} \lt= 20$:

$$
r_{lake(Secchi)} = r_{lake(Chla)}^{0.38}
$$

#### "Reference" conditions

$$
r_{lake(Secchi)} = r_{lake(Chla)}^{1.46}
$$










































