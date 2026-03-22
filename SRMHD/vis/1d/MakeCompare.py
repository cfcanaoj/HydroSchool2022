import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
import argparse


def read_textdata(file):
    with open(file, "r") as f:
        time = float(f.readline().split()[-1])

    arr = np.loadtxt(file)

    x = arr[:, 0]

    names = ["rho", "vx", "vy", "vz", "pre", "Bx", "By", "Bz", "sc"]
    data_dict = {name: arr[:, i+1] for i, name in enumerate(names)}

    return x, time, data_dict


#------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Create a comparison shapshot from multiple simulation outputs.",
    usage="python3 MakeCompare.py [step] [dir1] [dir2]",
    epilog=(
        "Example:\n"
        "  python3 MakeCompare.py 1 cw\n"
        "  python3 MakeCompare.py 1 hll hllc\n"
    ),
    formatter_class=argparse.RawTextHelpFormatter
)
parser.add_argument("step", type=int, help="step number")
parser.add_argument("dirnames", nargs="+", help="one or more directories containing snapshot files")
parser.add_argument("--anasol", default=None, help="exact solution file (optional)")
parser.add_argument("--xlim", nargs=2, type=float, default=[0.0, 1.0], metavar=("XMIN", "XMAX"))
args = parser.parse_args()
step       = args.step
dirnames   = args.dirnames

fig = plt.figure(figsize=(10, 10))
gs = gridspec.GridSpec(2, 2, wspace=0.30, hspace=0.2)
ax = []
for i in range(4):
    ax.append(plt.subplot(gs[i]))

ylabel = [r"$\rho$", r"$v_x$", r"$v_y$", r"$B_y$"]
for i, ax0 in enumerate(ax):
    ax0.minorticks_on()
    ax0.set_xlim(args.xlim[0], args.xlim[1])
    ax0.set_ylabel(ylabel[i])

ax[2].set_xlabel(r"$x$")
ax[3].set_xlabel(r"$x$")

#------------------------------------------------------------
# exact solution (optional)
#------------------------------------------------------------
has_anasol = False
if args.anasol is not None and os.path.exists(args.anasol):
    anasol = np.loadtxt(args.anasol)
    tout_ana = 0.1
    x_ana = anasol[:, 0]

    den_ana = np.concatenate(([anasol[0, 1]], anasol[:, 1], [anasol[-1, 1]]))
    vx_ana  = np.concatenate(([anasol[0, 2]], anasol[:, 2], [anasol[-1, 2]]))
    vy_ana  = np.concatenate(([anasol[0, 3]], anasol[:, 3], [anasol[-1, 3]]))
    by_ana  = np.concatenate(([anasol[0, 7]], anasol[:, 7], [anasol[-1, 7]]))
    has_anasol = True

#------------------------------------------------------------
# plot
#------------------------------------------------------------
time_ref = None
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

for idir, dirname in enumerate(dirnames):
    foutname = os.path.join(dirname, "snap%05d.dat" % step)
    print("reading " + foutname)

    x, time, data = read_textdata(foutname)

    if time_ref is None:
        time_ref = time

        if has_anasol:
            xlo, xhi = args.xlim
            x_ana1 = np.concatenate(([xlo], x_ana * time_ref / tout_ana, [xhi]))
            ax[0].plot(x_ana1, den_ana, '-', c="k", lw=1.5, label="exact")
            ax[1].plot(x_ana1, vx_ana,  '-', c="k", lw=1.5, label="exact")
            ax[2].plot(x_ana1, vy_ana,  '-', c="k", lw=1.5, label="exact")
            ax[3].plot(x_ana1, by_ana,  '-', c="k", lw=1.5, label="exact")

    color = colors[idir % len(colors)]
    label = os.path.basename(os.path.normpath(dirname))

    ax[0].plot(x, data['rho'], 'o-', c=color, mfc="none", label=label)
    ax[1].plot(x, data['vx'],  'o-', c=color, mfc="none", label=label)
    ax[2].plot(x, data['vy'],  'o-', c=color, mfc="none", label=label)
    ax[3].plot(x, data['By'],  'o-', c=color, mfc="none", label=label)

ax[0].text(
    0.5, 1.02, r"$\mathrm{time} = %.2f$" % (time_ref),
    horizontalalignment="center",
    transform=ax[0].transAxes
)

ax[0].legend(fontsize=9)

#------------------------------------------------------------
# output
#------------------------------------------------------------
if len(dirnames) > 1:
    outroot = "compare"
    for dirname in dirnames:
        outroot += "_" + os.path.basename(os.path.normpath(dirname))
else:
    outroot = os.path.basename(os.path.normpath(dirnames[0]))

for format_fig in ["pdf", "png"]:
    outdir = format_fig + "file"
    os.makedirs(outdir, exist_ok=True)

    outputfile = os.path.join(outdir, "snap%05d_%s.%s" % (step, outroot, format_fig))
    print("making plot file", outputfile)
    plt.savefig(outputfile, bbox_inches="tight")

plt.show()
# plt.close()
