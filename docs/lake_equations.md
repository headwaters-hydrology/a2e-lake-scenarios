## Lake equations for changes in input loads
These are the simplified lake equations solved for the change ratio as defined as the new input load (from reductions to to mitigations) to the original input load. These are derived from Abell et al 2019 and Abell et al 2020.

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

$r_{*}$ is the load change ratio, $z_{max}$ is the max depth of the lake in meters, and $\tau$ is the residence time of the lake in years.

Else when $\hspace{2mm} z_{max} \lt= 7.5$:

$$
r_{lake} = r_{in}
$$

#### "Reference" conditions
$$
\text{when} \hspace{2mm} z_{max} \gt 7.5: \\
r_{lake} = r_{in}^{1/b_{ref}}\\
\text{where...}\\
b_{ref} = 1 + 0.27 \tau^{0.29}\\
\text{else when} \hspace{2mm} z_{max} \lt= 7.5: \\
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
$$
r_{lake} = r_{in}^{1.25}\\
\text{Assuming a fixed ratio between } TN_{lake}/TP_{lake}
$$

#### "Reference" conditions
$$
r_{lake} = r_{in}^{1.24}\\
\text{Assuming a fixed ratio between } TN_{lake}/TP_{lake}
$$

### Secchi
#### Current conditions
$$
\text{when} \hspace{2mm} z_{max} \gt 20: \\
r_{lake} = r_{in}^{0.9}\\
\text{else when} \hspace{2mm} z_{max} \lt= 20: \\
r_{lake} = r_{in}^{0.38}
$$

#### "Reference" conditions
$$
r_{lake} = r_{in}^{1.46}
$$










































