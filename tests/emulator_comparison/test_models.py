import numpy as np
import matplotlib.pyplot as plt


# Load emulator data
yfull = np.loadtxt("./emulator/yFull.txt")
yshort = np.loadtxt("./emulator/yShort.txt")
logk = np.loadtxt("./emulator/logk.txt")
params_ce = np.loadtxt("./emulator/params_ce.txt")
k = 10**logk

# Load "ground truth" (high kmax, high k-per-log-interval)
rfs_target = [
    np.loadtxt(f"./output_kmax50_klogint1000/redTime_M{i+1:03d}.dat") for i in range(32)
]

# Load data to compare
rfs = [np.loadtxt(f"./output/redTime_M{i+1:03d}.dat") for i in range(32)]


def get_noneutrino_lin(k, rf, h):
    return rf[-128:, 3] / h**3 / (2 * np.pi**2) * k**1.5


def get_noneutrino_nlin(k, rf, h):
    return rf[-128:, 7] / h**3 / (2 * np.pi**2) * k**1.5


def get_neutrino_nlin(k, rf, h, om_nu, om_m):
    lin = get_noneutrino_lin(k, rf, h)
    nlin = get_noneutrino_nlin(k, rf, h)

    lin_nu = rf[-128:, 6] / h**3 / (2 * np.pi**2) * k**1.5

    # neutrino correction
    trans_p = np.sqrt(lin_nu / lin)
    beta_p = trans_p * (om_nu / om_m)
    f = 1.0 - om_nu / om_m + beta_p
    nlin = nlin * f**2
    return nlin


def test_neutrinoless():
    fig, axes = plt.subplots(2, 1, figsize=(7, 7))
    axes[0].set_title("z=0 emulator points vs redTime", y=0.95, va="top")
    axes[1].set_title("z=0 redTime vs redTime target", y=0.95, va="top")
    for ax in axes:
        ax.axhline(0, color="black")
        ax.axhspan(-0.001, 0.001, color="black", alpha=0.1)

    cmap = plt.matplotlib.colormaps["viridis"]

    for i in range(10):
        color = cmap(i / 9)
        h = params_ce.T[i][3]
        nlin_emulator = 10 ** yfull[-351:, i]

        # load data, interpolate to emulator k-values
        kcode = rfs[i][-128:, 0] * h  # correct to physical units
        nlin = get_noneutrino_nlin(kcode, rfs[i], h)
        nlin_int = 10 ** np.interp(logk[:40], np.log10(kcode), np.log10(nlin))

        # load truth data, interpolate to emulator k-values
        kcode_target = rfs_target[i][-128:, 0] * h  # correct to physical units
        nlin_target = get_noneutrino_nlin(kcode_target, rfs_target[i], h)
        nlin_target_int = 10 ** np.interp(
            logk[:40], np.log10(kcode_target), np.log10(nlin_target)
        )
        assert kcode.shape == kcode_target.shape

        axes[0].plot(
            k[:40],
            nlin_int / nlin_emulator[:40] - 1,
            color=color,
            linewidth=2,
        )
        axes[0].plot(
            k[:40],
            nlin_target_int / nlin_emulator[:40] - 1,
            color=color,
            linewidth=1,
            linestyle="dashed",
        )
        axes[1].plot(kcode, nlin / nlin_target - 1, color=color)

        kmask = kcode < 1e-1
        rel_diff = np.abs(nlin[kmask] / nlin_target[kmask] - 1)
        assert np.max(rel_diff) < 0.001
        assert np.all(np.isclose(kcode, kcode_target))

    axes[0].plot(
        [],
        color="black",
        linewidth=2,
        label="output",
    )
    axes[0].plot(
        [],
        color="black",
        linewidth=1,
        linestyle="dashed",
        label="targetoutput",
    )
    axes[0].legend(frameon=False)

    axes[0].set(xlabel="k (1/Mpc)", ylabel="code/emu - 1", xscale="log")
    axes[1].set(xlabel="k (1/Mpc)", ylabel="code/target - 1", xscale="log")
    fig.tight_layout()
    fig.savefig("neutrinoless_comparison.pdf", bbox_inches="tight")


def test_neutrinos():
    fig, axes = plt.subplots(2, 1, figsize=(7, 7))
    axes[0].set_title("z=0 emulator points vs inputs")
    axes[1].set_title("z=0 redTime vs redTime target", y=0.95, va="top")
    for ax in axes:
        ax.axhline(0, color="black")
        ax.axhspan(-0.001, 0.001, color="black", alpha=0.1)

    cmap = plt.matplotlib.colormaps["viridis"]

    for i in range(11, 32):
        color = cmap((i - 11) / (31 - 11))
        h = params_ce.T[i][3]
        om_nu = params_ce.T[i][7]
        om_m = params_ce.T[i][0]
        nlin_emulator = 10 ** yfull[-351:, i]

        # load data, interpolate to emulator k-values
        kcode = rfs[i][-128:, 0] * h  # correct to physical units
        nlin = get_neutrino_nlin(kcode, rfs[i], h, om_nu, om_m)
        nlin_int = 10 ** np.interp(logk[:40], np.log10(kcode), np.log10(nlin))

        # load truth data, interpolate to emulator k-values
        kcode_target = rfs_target[i][-128:, 0] * h  # correct to physical units
        nlin_target = get_neutrino_nlin(kcode_target, rfs_target[i], h, om_nu, om_m)
        nlin_target_int = 10 ** np.interp(
            logk[:40], np.log10(kcode_target), np.log10(nlin_target)
        )

        assert kcode.shape == kcode_target.shape
        assert np.all(np.isclose(kcode, kcode_target))

        axes[0].plot(
            k[:40], nlin_int / nlin_emulator[:40] - 1, color=color, linewidth=2
        )
        axes[0].plot(
            k[:40],
            nlin_target_int / nlin_emulator[:40] - 1,
            color=color,
            linewidth=1,
            linestyle="dashed",
        )
        axes[1].plot(kcode, nlin / nlin_target - 1, color=color)

        kmask = kcode < 1e-1
        rel_diff = np.abs(nlin[kmask] / nlin_target[kmask] - 1)
        assert np.max(rel_diff) < 0.005
        assert np.quantile(rel_diff, 0.95) < 0.001

    axes[0].plot(
        [],
        color="black",
        linewidth=2,
        label="output",
    )
    axes[0].plot(
        [],
        color="black",
        linewidth=1,
        linestyle="dashed",
        label="targetoutput",
    )
    axes[0].legend(frameon=False)

    axes[0].set(xlabel="k (1/Mpc)", ylabel="code/emu - 1", xscale="log")
    axes[1].set(xlabel="k (1/Mpc)", ylabel="code/target - 1", xscale="log")
    fig.tight_layout()
    fig.savefig("neutrino_comparison.pdf", bbox_inches="tight")
