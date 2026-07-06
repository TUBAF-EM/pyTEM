The package holds some functions for importing data for modelling IP effects, all based upon [pyGIMLi](https://pygimli.org).

Main entry, however, is the class TDIP, to be imported by

```python
from pyTEM import TEM
```

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




You