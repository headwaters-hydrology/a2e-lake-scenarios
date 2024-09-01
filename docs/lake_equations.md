## Lake equations for changes in input concentrations
These are the simplified lake equations solved for the change ratio as defined as the new input concentrations (from reductions to to mitigations) to the original input concentrations. These are derived from Abell et al 2019 and Abell et al 2020.

### TP
#### Current conditions
When $\hspace{2mm} z_{max} \gt 7.5$:

$$
r_{lake} = r_{in}^{1/b}
$$

where...

$$
b = 1 + 0.44 \tau^{0.13}
$$

$r_{*}$ is the concentration change ratio, $z_{max}$ is the max depth of the lake in meters, and $\tau$ is the residence time of the lake in years.

Else when $\hspace{2mm} z_{max} \lt= 7.5$:

$$
r_{lake} = r_{in}
$$

#### "Reference" conditions
When $\hspace{2mm} z_{max} \gt 7.5$:

$$
r_{lake} = r_{in}^{1/b_{ref}}
$$

where...

$$
b_{ref} = 1 + 0.27 \tau^{0.29}
$$

Else when $\hspace{2mm} z_{max} \lt= 7.5$:

$$
r_{lake} = r_{in}
$$

### TN
#### Current conditions

$$
r_{lake} = r_{in}^{0.54}
$$

#### "Reference" conditions

$$
r_{lake} = r_{in}^{0.81}
$$

### Chla
#### Current conditions
Assuming a fixed ratio between $TN_{lake}/TP_{lake}$:

$$
r_{lake} = r_{in}^{1.25}
$$

#### "Reference" conditions
Assuming a fixed ratio between $TN_{lake}/TP_{lake}$:

$$
r_{lake} = r_{in}^{1.24}
$$

### Secchi
#### Current conditions
When $\hspace{2mm} z_{max} \gt 20$:

$$
r_{lake} = r_{in}^{0.9}
$$

Else when $\hspace{2mm} z_{max} \lt= 20$:

$$
r_{lake} = r_{in}^{0.38}
$$

#### "Reference" conditions

$$
r_{lake} = r_{in}^{1.46}
$$










































