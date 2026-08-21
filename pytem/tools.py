"""Tools for TEM data import and processing."""
import numpy as np
import pandas as pd
from pathlib import Path

mu0 = 4e-7 * np.pi


def rhoa(t, dbzdt, m=1.0):
    return mu0 / np.pi * (mu0 * m / 20.0)**(2/3) * np.abs(dbzdt)**(-2/3) * t**(-5/3)

# dr = dr/du * du = -2/3 r / u *du => |dr/r| = 2/3 |dr/r|


def skinDepthTEM(t, rho):
    """Return diffusion (skin) depth of TEM."""
    return np.sqrt(2 * t * rho / mu0)


def bandpass(inp, p_dict):
    """Butterworth-type filter (implemented from simpegEM1D.Waveforms.py)."""
    cutofffreq = 1e8  # Determined empirically for TEM-FAST
    h = (1 + 1j*p_dict["freq"]/cutofffreq)**-1
    h *= (1 + 1j*p_dict["freq"]/3e5)**-1
    p_dict["EM"] *= h[:, None]

def readXYZfile(filename):
    """Read XYZ (workbench export) file.

    Parameters
    ----------
    filename : str
        Path to the XYZ file.

    Returns
    -------
    data : pandas.DataFrame
        DataFrame containing the XYZ data.
    header : dict
        Dictionary containing the header information.
    """
    with Path(filename).open() as fid:
        name = None
        header = {}
        lines = fid.readlines()
        for i, line in enumerate(lines):
            if line[0] != "/":
                break
            if name is None:
                name = line[1:].rstrip("\n")
            else:
                header[name] = line[1:]
                name = None

    nams = lines[i-1][1:].replace("[", "").replace("]", "").split()
    data = pd.read_csv(filename, delimiter=r"\s+", skiprows=i,
                       names=nams)
                    #    names=lines[i-1][1:].replace("[", "").split())
    return data, header
