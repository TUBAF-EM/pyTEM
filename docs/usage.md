# Usage

The package holds some functions for importing data for modelling IP effects, all based upon [pyGIMLi](https://pygimli.org).

Main entry, however, is the class TDIP, to be imported by

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

You can display the positions of the individual soundings by

```python
ax = p.showPositions()
```

The resulting axes object can be used to underlay any map.

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
