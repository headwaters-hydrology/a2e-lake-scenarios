# olw3-scenarios

### Processing procedure
1. Get LAWA river WQ sample site data: lawa-wq_data.
2. Delineate catchments above all LAWA sites.
3. Process reference conc at all LAWA sites from the data from Rich.
4. Process the land cover from the LCDB and the special Dairy and Sheep/Beef layers/data from Rich and Ton.
   * Yield values are taken from a combo of MS and Rich and reductions are from Rich.
   * All values are assigned land covers at their original classification, then aggregated to the more generic land cover classes for this project. This mean that each catchment will have slightly different yield values for the generic classes.
   * Ultimately, the relative proportions of each class within each catchment is needed for the final web app calcs.


### E.coli land area proportioning equations
This describes some mass balance equations regarding areas of a catchment with two land uses (sheep and beef (S) and dairy (D)) with two drainage classes. This also includes a state change from the current state (1) to the user-modified (scenario) state (2). There are several more land uses and a set of elevation classes, but this is meant as a simple generalised example that can be applied to the full set of classes.

The user can only change the areas associated with the land uses and not the subclasses of elevation and drainage. We need to know the areas for all the subclasses of all the land uses, so we need to make some assumptions about how to distribute the areas within the subclasses when the user changes the land use areas.

For the following examples, we will assume that there are only four land uses as this will make the examples simpler, but this method will work for any number of land uses. These land uses include Dairy ($D$), Sheep and beef ($S$), Exotic forest/Forestry ($F$), and Native vegetation ($N$). There are two drainage classes (i.e. $well$ and $poor$) and two elevation classes (i.e. $high$ and $low$).

The concept is that any land uses (and associated subclasses) that have decreased in size should be added up and those land uses that have increased in size should then receive those added up areas. Since there is no clear way to know what specific land uses are converted to other land uses, the new land areas must be proportioned by newly added land area.

Some examples of the basic mass balance equations are shown below. The superscript reflects the land use type and the subscript reflects the subclass. The term $catch$ refers to the entire catchment.

$$
% \tag{1}
A_{catch}^{catch} = A_{well,high}^{catch} + A_{well,low}^{catch} + A_{poor,high}^{catch} + A_{poor,low}^{catch}
$$

$$
% \tag{2}
1 = \frac{A_{well,high}^{catch}}{A_{catch}^{catch}} + \frac{A_{well,low}^{catch}}{A_{catch}^{catch}} + \frac{A_{poor,high}^{catch}}{A_{catch}^{catch}} + \frac{A_{poor,low}^{catch}}{A_{catch}^{catch}}
$$

Each of the $A_{*}^{catch}$ will be associated with all land uses. For example, if we only have two land uses (Dairy (D) and Sheep and beef (S)), then the generalised equation would be:

$$
% \tag{3}
A_{*}^{catch} = A_{*}^{D} + A_{*}^{S} + A_{*}^{F} + A_{*}^{N}
$$

What we know are the ratios in state 1 (current state) for equations (2) and (3), but we don't have this information in state 2 (scenario state). In state 2, we know the aggregate ratio of the new area in each land use. For example, for Dairy we would know $R_{catch}^{D}$:

$$
% \tag{4}
R_{catch(2)}^{D} = \frac{A_{well,high(2)}^{D} + A_{well,low(2)}^{D} + A_{poor,high(2)}^{D} + A_{poor,low(2)}^{D}}{A_{catch}^{catch}}
$$

But we do not know the individual ratios (e.g. $A_{well,high(2)}^{D}$).

The calculation process initially involves estimating the areas that have decreased in size from state 1 to 2. When $R_{catch(2)}^{*} < R_{catch(1)}^{*}$, subtract the ratios and multiply by all areas associated with that land use. For example, if $R_{catch(1)}^{D}$ was 0.35 and $R_{catch(2)}^{D}$ is 0.2, then the lost areas for all the classes in the dairy land use would be:

$$
% \tag{5}
\Delta A_{well,high}^{D} = (R_{catch(1)}^{D} - R_{catch(2)}^{D})A_{well,high(1)}^{D}
$$

$$
% \tag{6}
\Delta A_{well,low}^{D} = (R_{catch(1)}^{D} - R_{catch(2)}^{D})A_{well,low(1)}^{D}
$$

$$
% \tag{7}
\Delta A_{poor,high}^{D} = (R_{catch(1)}^{D} - R_{catch(2)}^{D})A_{poor,high(1)}^{D}
$$

$$
% \tag{8}
\Delta A_{poor,low}^{D} = (R_{catch(1)}^{D} - R_{catch(2)}^{D})A_{well,poor(1)}^{D}
$$

and $(R_{catch(1)}^{D} - R_{catch(2)}^{D})$ would be 0.15.

This would be performed for every land use that has decreased in size. Then those four subclasses (i.e. $(well,high)$, $(well,low)$, $(poor,high)$, and $(poor,low)$) would have their areas aggregated. For example, if both D and S lost areas, then:

$$
% \tag{9}
\Delta A_{well,high}^{catch} = \Delta A_{well,high}^{D} + \Delta A_{well,high}^{S}
$$

$$
% \tag{10}
\Delta A_{well,low}^{catch} = \Delta A_{well,low}^{D} + \Delta A_{well,low}^{S}
$$

$$
% \tag{11}
\Delta A_{poor,high}^{catch} = \Delta A_{poor,high}^{D} + \Delta A_{poor,high}^{S}
$$

$$
% \tag{12}
\Delta A_{poor,low}^{catch} = \Delta A_{poor,low}^{D} + \Delta A_{poor,low}^{S}
$$

This gives us the four subclasses with the lost areas. These lost areas need to be added to the land uses that have increased from state 1 to 2. Similar to the land uses that lost areas, when $R_{catch(2)}^{*} > R_{catch(1)}^{*}$, subtract the ratios and multiply by all areas associated with that land use. In our above example where $(R_{catch(1)}^{D} - R_{catch(2)}^{D})$ is 0.15 and $(R_{catch(1)}^{S} - R_{catch(2)}^{S})$ is 0.05 (a total change of 0.2), if Exotic Forest ($R_{catch}^{F}$) increased from 0.1 to 0.25 and native vegetation ($R_{catch}^{N}$) increased from 0.1 to 0.15, then the extra areas that would need to be added to the two land uses would be:

$$
% \tag{13}
A_{well,high(2)}^{F} = \frac{(R_{catch(2)}^{F} - R_{catch(1)}^{F})}{(R_{catch(1)}^{D} - R_{catch(2)}^{D}) + (R_{catch(1)}^{S} - R_{catch(2)}^{S})} \Delta A_{well,high}^{catch} A_{well,high(1)}^{F}
$$

With our actual numbers we get:

$$
% \tag{14}
A_{well,high(2)}^{F} = \frac{(0.25 - 0.1)}{(0.35 - 0.2) + (0.45 - 0.4)} \Delta A_{well,high}^{catch} A_{well,high(1)}^{F}
$$

This would be applied to all land uses and subclasses that have increased in area. Consequently, we would know the areas for all the land uses and subclasses for state 1 and 2. This would allow us to calculate the change in yields from state 1 to state 2 using Ton's model.


