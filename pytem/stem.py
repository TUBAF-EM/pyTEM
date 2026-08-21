"""Tools for reading and processing sTEM data."""
from pathlib import Path
import numpy as np


def readGEXFile(fname="sTEM.gex"):
    """Read GEX file (a rather general function)."""
    with Path(fname).open() as fid:
        lines = fid.readlines()
        for i, line in enumerate(lines):
            lines[i] = line.replace("=", "= ")

    out = {}
    putin = out
    for line in lines:
        line = line.replace("\n", "").replace("\r", "")
        if len(line) > 2 and line[0] == "[" and line[-1] == "]":
            sect = line[1:-1]
            if sect != "General":
                putin = out[sect] = {}

        if "=" in line:
            nam, val = line.split("=")
            try:
                putin[nam] = np.fromstring(val, dtype=float, sep=" ")
            except ValueError:
                putin[nam] = val

    return out

def collectNumData(dic, name, start=1, stop=100, num=0):
    """Collect data from dictionary stored in numerically named keys."""
    if num == 0:
        for num in range(1, 4):
            if f"{name}{start:0{num}d}" in dic:
                break
    i = start
    col = []
    for i in range(start, stop):
        key = f"{name}{i:0{num}d}"
        if key in dic:
            col.append(dic[key])
        else:
            break

    return np.array(col)


def readSettings(filename="sTEM.gex", ramptime=None):
    """Read settings."""
    out = readGEXFile(filename)
    cfg = {}
    cfg["rxpos"] = out["RxCoilPosition1"]
    cfg["txpos"] = out["TxCoilPosition1"]
    if "TxLoopPoint1" in out:
        txp = collectNumData(out, "TxLoopPoint")
        cfg["tx"], cfg["ty"] = txp[:, 0], txp[:, 1]
        xx = np.hstack([cfg["tx"], cfg["tx"][0]])
        yy = np.hstack([cfg["ty"], cfg["ty"][0]])
        cfg.setdefault("txarea", np.abs(sum(xx[:-1]*yy[1:]) -
                                        sum(xx[1:]*yy[:-1]))/2)
    else:
        dia = out["TxLoopDiameter"]
        cfg.setdefault("txarea", dia**2*np.pi/4)
        nL = 8
        cfg["tx"] = np.sin(np.arange(nL) * 2 * np.pi / nL) * dia / 2
        cfg["ty"] = np.cos(np.arange(nL) * 2 * np.pi / nL) * dia / 2

    if "WaveformLMPoint01" in out: # dual mode
        bla = collectNumData(out, "WaveformLMPoint", num=2)
        cfg["tL"], cfg["vL"] = bla[:, 0], bla[:, 1]
        bla = collectNumData(out, "WaveformHMPoint", num=2)
        cfg["tH"], cfg["vH"] = bla[:, 0], bla[:, 1]
        cfg["timeL"] = collectNumData(out, "GateTimeLM")[:, 0]
        cfg["timeH"] = collectNumData(out, "GateTimeHM")[:, 0]
    else:
        if "WaveformPoint" in out: # single mode
            bla = collectNumData(out, "WaveformPoint", num=2)
            cfg["t"], cfg["v"] = bla[:, 0], bla[:, 1]
        else:
            if ramptime is None:
                print("No waveform found. Guessing.")
                ramptime = 5e-3

            cfg["t"] = np.array([-1e-4, 0, ramptime])
            cfg["v"] = np.array([1.0, 1.0, 0.0])

        cfg["time"] = collectNumData(out, "GateTime")[:, 0]

    return cfg

def readSettings1(filename="sTEM.gex"):
    with Path(filename).open() as fid:
        lines = fid.readlines()

    cfg = {}
    cfg["tL"], cfg["vL"] = np.genfromtxt(lines[20:60], usecols=[1, 2], unpack=True)
    cfg["tH"], cfg["vH"] = np.genfromtxt(lines[61:109], usecols=[1, 2], unpack=True)
    cfg["timeL"] = np.genfromtxt(lines[111:117], usecols=[1])
    cfg["timeH"] = np.genfromtxt(lines[119:141], usecols=[1])
    cfg["tx"], cfg["ty"] = np.genfromtxt(lines[15:19], usecols=[1, 2], unpack=True)
    cfg["rxpos"] = np.genfromtxt(lines[7:8], usecols=[1, 2, 3])
    cfg["txarea"] = np.genfromtxt(lines[14:15], usecols=[1])
    return cfg


def bandpass(inp, p_dict):
    """Butterworth-type filter (implemented from simpegEM1D.Waveforms.py)."""
    cutofffreq = 1e8  # Determined empirically for TEM-FAST
    h = (1 + 1j*p_dict["freq"]/cutofffreq)**-1
    h *= (1 + 1j*p_dict["freq"]/3e5)**-1
    p_dict["EM"] *= h[:, None]


# %%
if __name__ == "__main__":
    pass