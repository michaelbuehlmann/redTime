# imports
import numpy as np

params_list = ["om_m", "om_b", "s_8", "h", "n_s", "w_0", "-(w0+wa)^1/4", "om_nu"]
ranges_lower = np.array([0.12, 0.0215, 0.7, 0.55, 0.85, -1.3, 0.3, 0.0])
ranges_upper = np.array([0.155, 0.0235, 0.9, 0.85, 1.05, -0.7, 1.29, 0.01])
diff_ranges = ranges_upper - ranges_lower

# I/O
lhc = np.loadtxt("s-lhs.100.8_1")

with open("params_lhc.dat", "w") as f:
    # write preamble
    f.write("# Cosmological models (1 per line)\n")
    f.write("#\n")
    f.write("# Columns\n")
    f.write(
        "#model  omega_m omega_b s8       h       ns      w0       wa       omega_nu\n"
    )
    f.write("#\n")

    # write models
    for l in range(len(lhc)):
        vals = lhc[l] * diff_ranges + ranges_lower
        f.write("M{:03d}".format(l + 1))
        for j in range(8):
            if j == 6:
                val_wa = -((vals[6]) ** 4) - vals[5]
                f.write("  " + str(val_wa))
            else:
                f.write("  " + str(vals[j]))
        f.write("\n")
