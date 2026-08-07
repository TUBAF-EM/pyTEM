"""TEM modeling classes."""
import numpy as np
from empymod import bipole
import pygimli as pg
from .tools import bandpass
from .stem import readSettings

class TEMBlockModelling(pg.frameworks.Block1DModelling):
    """TEM modelling class using block 1D discretization."""

    def __init__(self, **kwargs):
        """Initialize instance.

        Parameters
        ----------
        cfg : str|dict
            configuration file or dictionary, containing
            - t, v : signal waveform
            - time : time (gate midpoint) vector
            - tx, ty : transmitter polygon
            - rxpos : receiver position
            - txarea : transmitter area
        """
        cfg = kwargs.get("cfg", {})
        super().__init__(**kwargs)
        if isinstance(cfg, str):
            cfg = readSettings(cfg)

        if "t" in cfg and "v" in cfg:
            self.signal = {'nodes': cfg["t"], 'amplitudes': cfg["v"], 'signal': 1}
        else:
            self.signal = -1

        self.time = cfg["time"]
        # TODO: compute txarea from points if not present

        self.kw = dict(
            src=[[cfg["tx"][-1], *cfg["tx"]], # x1
                 [*cfg["tx"], cfg["tx"][0]], # x2
                 [cfg["ty"][-1], *cfg["ty"]], # y1
                 [*cfg["ty"], cfg["ty"][0]], # y2
            0, 0],
            strength=1/cfg["txarea"],
            verb=0,
            rec = np.concatenate([cfg["rxpos"], [0, 90]]),       # Receiver at the origin, vertical.
            mrec="b",                   # Receiver: dB/dt
            srcpts=3,                   # Approx. the finite dip. with 3 points.
            ftarg={"dlf": "key_81_2009"},  # Shorter, faster filters.
            htarg={"dlf": "key_101_2009", "pts_per_dec": -1},
            bandpass={"func": bandpass}
            )

    @property
    def t(self):
        """Return time vector."""
        return self.time

    def response(self, model):
        """Return model response."""
        thk = model[:self.nLayers-1]
        res = model[self.nLayers-1:]
        return bipole(
            depth=np.concatenate([[0], np.cumsum(np.atleast_1d(thk))]), # Depth-model.
            res=np.concatenate([[2e14], np.atleast_1d(res)]),      # Resistivity model.
            signal=self.signalL,
            freqtime=self.timeL,      # Wanted times.
            **self.kw).sum(axis=1)


class TEMBlockModellingDualMode(TEMBlockModelling):
    """TEM modelling class using block 1D discretization."""

    def __init__(self, **kwargs):
        """Initialize instance.

        Parameters
        ----------
        cfg : str|dict
            configuration file or dictionary, containing
            - tL, vL : low moment signal waveform
            - tH, vH : high moment signal waveform
            - timeL, timeH : time (gate midpoint) vectors
            - tx, ty : transmitter polygon
            - rxpos : receiver position
            - txarea : transmitter area
        """
        cfg = kwargs.pop("cfg", {})
        super().__init__(**kwargs)
        if isinstance(cfg, str):
            cfg = readSettings(cfg)

        self.signalL = {'nodes': cfg["tL"], 'amplitudes': cfg["vL"], 'signal': 1}
        self.signalH = {'nodes': cfg["tH"], 'amplitudes': cfg["vH"], 'signal': 1}
        self.timeL = cfg["timeL"]
        self.timeH = cfg["timeH"]


    def response(self, model):
        """Return model response."""
        thk = model[:self.nLayers-1]
        res = model[self.nLayers-1:]
        outL = bipole(
            depth=np.concatenate([[0], np.cumsum(np.atleast_1d(thk))]), # Depth-model.
            res=np.concatenate([[2e14], np.atleast_1d(res)]),      # Resistivity model.
            signal=self.signalL,
            freqtime=self.timeL,      # Wanted times.
            **self.kw)
        outH = bipole(
            depth=np.concatenate([[0], np.cumsum(np.atleast_1d(thk))]), # Depth-model.
            res=np.concatenate([[2e14], np.atleast_1d(res)]),      # Resistivity model.
            signal=self.signalH,
            freqtime=self.timeH,      # Wanted times.
            **self.kw)
        return np.concatenate([outL.sum(axis=1), outH.sum(axis=1)])


class TEMRhoModelling(pg.frameworks.MeshModelling):
    """TEM modelling class using fixed (smooth) 1D discretization."""

    def __init__(self, thk, **kwargs):
        """Initialize class instance.

        Parameters
        ----------
        cfg : str|dict
            configuration file or dictionary, containing
            - t, v : signal waveform
            - time : time (gate midpoint) vector
            - tx, ty : transmitter polygon
            - rxpos : receiver position
            - txarea : transmitter area
        """
        self.thk = thk
        cfg = kwargs.pop("cfg", "sTEM.gex")
        self.mesh_ = pg.meshtools.createMesh1D(len(thk)+1)
        super().__init__(mesh=self.mesh_)
        if isinstance(cfg, str):
            cfg = readSettings(cfg)

        if "t" in cfg and "v" in cfg:
            self.signal = {'nodes': cfg["t"], 'amplitudes': cfg["v"], 'signal': 1}
        else:
            self.signal = -1

        self.kw = dict(
            src=[[cfg["tx"][-1], *cfg["tx"]], # x1
                 [*cfg["tx"], cfg["tx"][0]], # x2
                 [cfg["ty"][-1], *cfg["ty"]], # y1
                 [*cfg["ty"], cfg["ty"][0]], # y2
                 0, 0],
            strength=1/cfg["txarea"],
            verb=0,
            depth = np.concatenate([[0], np.cumsum(np.atleast_1d(thk))]),
            rec = np.concatenate([cfg["rxpos"], [0, 90]]), # Receiver at the origin, vertical.
            mrec="b",                   # Receiver: dB/dt
            srcpts=3,                   # Approx. the finite dip. with 3 points.
            ftarg={"dlf": "key_81_2009"},  # Shorter, faster filters.
            bandpass={"func": bandpass}
            )
        if kwargs.pop("ht", False):
            self.kw["htarg"] = {"dlf": "key_101_2009", "pts_per_dec": -1}

    @property
    def t(self):
        """Return time vector."""
        return np.concatenate([self.timeL, self.timeH])

    def response(self, model):
        """Return model response."""
        return bipole(res=np.concatenate([[2e14], model]),
                   signal=self.signalL,
                   freqtime=self.timeL,      # Wanted times.
                   **self.kw).sum(axis=1)

    def createStartVector(self, data):
        """Create a starting vector."""
        return pg.Vector(len(self.thk)+1, 250.)

class TEMRhoModellingDualMode(TEMRhoModelling):
    """TEM modelling class using fixed (smooth) 1D discretization."""

    def __init__(self, thk, **kwargs):
        """Initialize class instance.

        Parameters
        ----------
        cfg : str|dict
            configuration file or dictionary, containing
            - tL, vL : low moment signal waveform
            - tH, vH : high moment signal waveform
            - timeL, timeH : time (gate midpoint) vectors
            - tx, ty : transmitter polygon
            - rxpos : receiver position
            - txarea : transmitter area
        """
        cfg = kwargs.get("cfg", "sTEM.gex")
        super().__init__(thk, **kwargs)
        if isinstance(cfg, str):
            cfg = readSettings(cfg)

        self.signalL = {'nodes': cfg["tL"], 'amplitudes': cfg["vL"], 'signal': 1}
        self.signalH = {'nodes': cfg["tH"], 'amplitudes': cfg["vH"], 'signal': 1}
        self.timeL = cfg["timeL"]
        self.timeH = cfg["timeH"]

    @property
    def t(self):
        """Return time vector."""
        return np.concatenate([self.timeL, self.timeH])

    def response(self, model):
        """Return model response."""
        return np.concatenate([
            bipole(res=np.concatenate([[2e14], model]),
                   signal=self.signalL,
                   freqtime=self.timeL,      # Wanted times.
                   **self.kw).sum(axis=1),
            bipole(res=np.concatenate([[2e14], model]),
                   signal=self.signalH,
                   freqtime=self.timeH,      # Wanted times.
                   **self.kw).sum(axis=1)])

# %%
if __name__ == "__main__":
    f = TEMRhoModelling(thk=np.arange(2, 28, 2), cfg="sTEM.gex")
    rho = pg.Vector(len(f.thk)+1, 100.)
    print(f(rho))
