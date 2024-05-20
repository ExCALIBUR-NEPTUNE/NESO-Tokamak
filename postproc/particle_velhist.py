import functools
import h5py
import numpy as np
import os.path
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from scipy.optimize import curve_fit

from common import get_run_dir, get_postproc_dir


def vhist_movie(run_lbl, fname="particle_trajectory.h5part", ylim=0.6, **hist_settings):
    # Default histogram settings
    nbins = hist_settings.pop("bins", 20)
    hrange = hist_settings.pop("range", (0, 7.5))
    hist_settings["lw"] = hist_settings.get("lw", 1)
    hist_settings["ec"] = hist_settings.get("ec", "grey")
    hist_settings["fc"] = hist_settings.get("fc", "red")
    hist_settings["alpha"] = hist_settings.get("alpha", 0.7)

    # Read data, sort into chronological order
    run_dir = get_run_dir(run_lbl)
    data = h5py.File(os.path.join(run_dir, fname), "r")
    keys = sorted(data.keys(), key=lambda x: int(x.split("#")[-1]))
    hist_norm = data[keys[0]]["VELOCITY_0"].size * (hrange[1] - hrange[0]) / nbins

    # Function to get masses and velocity magnitudes at a particular step
    def get_data(step):
        step_data = data[keys[step]]
        vmags = np.sqrt(
            step_data["VELOCITY_0"][()] * step_data["VELOCITY_0"][()]
            + step_data["VELOCITY_1"][()] * step_data["VELOCITY_1"][()]
            + step_data["VELOCITY_2"][()] * step_data["VELOCITY_2"][()]
        )
        # Assume mass=1 for all particles for now
        # masses = step_data["MASS"][()]
        masses = 1
        return masses, vmags

    # Set up axes
    fig, ax = plt.subplots()
    ax.set_ylim(top=ylim)
    ax.set_xlabel("|v|")
    ax.set_ylabel("Prob")

    # Init histogram
    HIST_BINS = np.linspace(hrange[0], hrange[1], nbins)
    _, _, bar_container = ax.hist(get_data(0)[1], HIST_BINS, **hist_settings)
    HIST_BIN_CENTERS = HIST_BINS[:-1] + np.diff(HIST_BINS) / 2

    # Init Maxwellian with fitted kT
    nbins_MB = 200
    MB_BINS = np.linspace(hrange[0], hrange[1], nbins_MB)
    mb_fitted_settings = dict(linestyle="-", color="blue", label="not set")
    (mb_fitted,) = ax.plot(MB_BINS, nbins_MB * [1.0], **mb_fitted_settings)
    #   Initial guess for fitted kT
    init_params = [0.1]
    #   Require kT > 0.1 (curve_fit doesn't like smaller values...)
    kT_bounds = ((0.1), (np.inf))
    legend = ax.legend()

    # Function to compute maxwellian. Mass kwarg can be scalar or array.
    def maxwellian(v, kT, m=1):
        return m * v / kT * np.exp(-m * np.square(v) / 2.0 / kT)

    # Function to update plot for each frame
    def animate(frame, bar_container):
        masses, vmags = get_data(frame)

        # Update histogram bar heights
        n, _ = np.histogram(vmags, HIST_BINS)
        pdf = n / hist_norm
        for count, rect in zip(pdf, bar_container.patches):
            rect.set_height(count)

        # Update 2D Maxwellian with kT fitted to the histogram
        maxwellian_passm = functools.partial(maxwellian, m=masses)
        fitted_params, _ = curve_fit(
            maxwellian_passm, HIST_BIN_CENTERS, pdf, p0=init_params, bounds=kT_bounds
        )
        # Update fitted kT val in legend
        mb_fitted.set_ydata(maxwellian_passm(MB_BINS, *fitted_params))
        legend.get_texts()[0].set_text(f"Maxwellian with kT = {fitted_params[0]:.2f}")

        return bar_container.patches

    # Generate and save the animation
    func = functools.partial(animate, bar_container=bar_container)
    anim = FuncAnimation(fig, func, len(keys), repeat=False, blit=True)
    anim_writer = FFMpegWriter(fps=15)
    anim.save(
        filename=os.path.join(get_postproc_dir(run_lbl), "vhist.mp4"),
        writer=anim_writer,
    )


# ==============================================================================

run_lbl = "implicit2D.prev"
hist_settings = dict()

vhist_movie(run_lbl, **hist_settings)
