import numpy as np
from pytem import TEM

tem = TEM("sTEM.xyz")
print(tem)

di = tem.f.kw.copy()
for i, li in enumerate(di["src"]):
    if hasattr(li, '__iter__'):
        for j, po in enumerate(li):
            li[j] = float(po)

    tem.f.kw["src"][i] = li

di["depth"] = [float(ff) for ff in di["depth"]]
di["rec"] = [float(ff) for ff in di["rec"]]
del di["bandpass"]
print(di)

import json
with open("sTEM-40x40m.json", "w") as f:
    json.dump(di, f, indent=4)# separators=("\n", ":"))