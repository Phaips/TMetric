import warnings
warnings.filterwarnings('ignore', category=FutureWarning)
import numpy as np
import matplotlib.pyplot as plt
import starfile
from scipy.spatial import cKDTree
from scipy.stats import gaussian_kde
ordinal = lambda n: f"{n}{'th' if 11<=n%100<=13 else {1:'st',2:'nd',3:'rd'}.get(n%10,'th')}"

star_path   = "path/to/particles.star"

k_max       = 8                        # number of NN
xlim        = (100, 300)                # adjust to your likeing
rdf_r_max, rdf_dr = 800, 20            # up to distance in Å, increment in Å
radii_nm = [15, 20, 25, 30, 35, 40]    # radii for neighbor‐counts plot in nm
radii_A  = [r * 10 for r in radii_nm]  # (nm → Å)

df = starfile.read(star_path)
coords = {
    t: g[["rlnCenteredCoordinateXAngst",
           "rlnCenteredCoordinateYAngst",
           "rlnCenteredCoordinateZAngst"]].values
    for t, g in df.groupby("rlnTomoName")
}
dist = {t: cKDTree(c).query(c, k_max+1)[0] for t, c in coords.items()}

neighbors = range(2, k_max+1)
cmap      = plt.get_cmap('tab10')
colors    = [cmap(i) for i in range(len(neighbors))]

for tomo, dmat in dist.items():
    coords_t = coords[tomo]
    tree     = cKDTree(coords_t)
    fig, ax  = plt.subplots(1, 4, figsize=(20, 4))
    labels0 = []
    for j, col in zip(neighbors, colors):
        vals = dmat[:, j]
        n, bins = np.histogram(vals, bins=100)
        centers = 0.5 * (bins[1:] + bins[:-1])
        peak0 = centers[n.argmax()]
        std   = vals.std()
        ax[0].hist(vals, bins=100, alpha=0.5, edgecolor='k', color=col)
        labels0.append(f"{ordinal(j-1)} NN: {peak0:.1f}±{std:.1f}Å")
    ax[0].set(title=tomo, xlabel="Distance (Å)", ylabel="Count", xlim=xlim)
    ax[0].legend(labels0, loc='upper right')
    labels1 = []
    x_vals = np.linspace(xlim[0], xlim[1], 200)
    for j, col in zip(neighbors, colors):
        vals = dmat[:, j]
        y    = gaussian_kde(vals)(x_vals)
        peak1 = x_vals[y.argmax()]
        std   = vals.std()
        ax[1].hist(vals, bins=100, density=True, alpha=0.3, edgecolor='k', color=col)
        ax[1].plot(x_vals, y, color=col)
        labels1.append(f"{ordinal(j-1)} NN: {peak1:.1f}±{std:.1f}Å")
    ax[1].set(title=tomo, xlabel="Distance (Å)", ylabel="Density", xlim=xlim)
    ax[1].legend(labels1, loc='upper right')
    counts = []
    for r_A in radii_A:
        nb_lists = tree.query_ball_point(coords_t, r_A)
        counts.append([len(nb)-1 for nb in nb_lists])
    positions = np.arange(len(radii_nm)) * 2
    bp = ax[2].boxplot(
        counts,
        positions=positions,
        widths=0.6,
        patch_artist=True,
        showfliers=False,
        manage_ticks=False
    )
    for patch in bp['boxes']:
        patch.set_facecolor('gray')
        patch.set_edgecolor('k')
        patch.set_alpha(0.8)
    ax[2].set_xticks(positions)
    ax[2].set_xticklabels([f"{r} nm" for r in radii_nm], fontweight='bold')
    ax[2].set(title="Neighbor Counts", xlabel="Radius", ylabel="Count", ylim=(0,140))

    bins = np.arange(0, rdf_r_max + rdf_dr, rdf_dr)
    hist = np.zeros(len(bins)-1)
    for c0 in coords_t:
        idx   = tree.query_ball_point(c0, rdf_r_max)
        dists = np.linalg.norm(coords_t[idx] - c0, axis=1)
        hist += np.histogram(dists[dists>0], bins=bins)[0]
    shellA = (4/3) * np.pi * (bins[1:]**3 - bins[:-1]**3)
    local_density = hist / (2 * len(coords_t) * (shellA / 1e7))
    r_centers     = 0.5 * (bins[:-1] + bins[1:])
    ax[3].plot(r_centers, local_density)
    ax[3].set(title="RDF", xlabel="Distance (Å)", ylabel="Local density")

    plt.tight_layout()
    plt.show()
