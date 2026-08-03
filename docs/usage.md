# Usage

The package holds some functions for importing data for modelling IP effects, all based upon [pyGIMLi](https://pygimli.org).

Main entry, however, is the class TEM, to be imported by

```python
from pyTEM import TEM
```

## Loading data

It can be initialized empty or directly with a data file, up to now,
only xyz (with accompagnying gex) and usf files are supported.

```python
p = TEM("mydata.xyz")
```

For xyz files, it tries to find a `.gex` file with the same name.
You can also specify this by

```python
p = TEM("mydata.xyz", cfg="myconf.gex")
```

By using

```python
print(p)
```

As a result, you will see something like

```
TEM profile with 20 soundings
47 times (1.245463e-05-0.01112942)
```

Typing `p.cfg` will show you a dictionary with all measuring parameters like

* transmitter position and area
* receiver positions
* time gates
* transmitter waveforms


## Transmitter waveform

The transmitter waveform(s) can be shown by

```python
ax = p.showWaveform();
```

Usually you mainly see the tx repetition rate and the on-switch in µs.
To have a closer look at the off-switch ramp, you need to adjust the axis a bit, e.g.

```python
ax.set_xlim(0, 20)
```

## Data

### Positions

You can display the positions of the individual soundings by

```python
ax = p.showPositions()
```

Any `show` commant returns a [matplotlib](https://matplotlib.org) axes object that can be modified freely (limits, labels) or be used to underlay any map, e.g. `underlayBKGMap` from pyGIMLi.

You may want to have a spatial look at some data, e.g. by plotting the apparent resistivity of a distinct time window (number `nt`)

```python
ax = tem.showRhoa(nt=10) # automatic
ax = tem.showRhoa(nt=10, cMin=20, cMax=500, orientation="vertical") # colorbar
```

If you do not specify `nt`, it will plot all apparent resistivities next to each other (x axis) and the time gates on the y axis.

```python
ax = tem.showRhoa() # cMin/cMax etc.
```

This helps to get an overview on the general variability and noise level, but can only be interpreted geometrically if the data are acquired along a profile.

### Soundings

You may have a look at single soundings as transient

```python
p.showSounding(n=11)
```

or apparent resistivity

```python
p.showSounding(n=11, rhoa=True)
```

To show a number of soundings, you can use

```python
p.showSoundings()  # all
ax = p.showSoundings([0, 1, 5, 6]) # only a few
ax.legend()
```

### Filtering

As a result, you will see when the transients run into noise. You can filter the time series with

```python
p.filter(tmax=1e-3)
```

You can also filter some soundings by their number

```python
tem.filter(n=[15, 25, 26])
```

or extract some of them (a profile or a part of a profile)

```python
tem.filter(nmin=10, nmax=22)
```

Note that this number is persistent after reading (like it is shown by `showPositions()`) and does not change upon removing soundings.

## Inversion

A single sounding can be inverted by

```python
tem.invertSounding(20);
```

It will directly show the model and its response on top of the data in a two-column figure (two ax returned).

Up to now, a multi-layer (Occam) inversion is done.
You can specify the (fixed) thickness layer vector in the inversion (only for this call) call or directly upon the initialization (for all inversions) using `thk=`.

You can also invert all data using

```python
tem.invertAll();
```

Several options of constraining them towards each other are currently implemented, by default the result of one sounding is taken as starting model for the next.
The result is by default shown as 2D profile, but can also be shown by

```python
p.showResults()
```

Here, you can specify the colorbar limits (`cMin`/`cMax`) the colormap (`cMap`), log behaviour for the color (`logScale`) or the depth (`zLog`) and depth limits (`zMin`/`zMax`).

### Visualization in 3D

Export to VTK.

## Example notebooks

In the [github data folder](https://github.com/TUBAF-EM/pyTEM/tree/main/data) you will find some example data with notebooks demonstrating how to work with it