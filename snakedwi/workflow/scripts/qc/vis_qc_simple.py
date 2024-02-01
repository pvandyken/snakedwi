import matplotlib
import matplotlib.pyplot as plt
import nibabel as nib
import skimage as skimg
from nilearn import plotting

matplotlib.use("Agg")


with plt.ioff():
    background = nib.load(snakemake.input["img"])

    # Static Edge overlay
    ax_size = 2
    n_slices = 7
    fig = plt.figure(
        figsize=(ax_size * n_slices, ax_size * 4),
        layout="constrained",
        facecolor="k",
        edgecolor="k",
    )
    ax = fig.add_subplot(111)
    display = plotting.plot_img(
        background.slicer[..., 0],
        n_slices,
        display_mode="mosaic",
        black_bg=True,
        colorbar=False,
        annotate=False,
        cmap="gray",
        axes=ax,
    )
    ax.axis("off")
    fig.savefig(snakemake.output.png, dpi=600)
    plt.close(fig)

