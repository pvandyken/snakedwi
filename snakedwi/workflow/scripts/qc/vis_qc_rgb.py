from nilearn import plotting
from matplotlib import pyplot as plt
import io
import nibabel as nib
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

dpi = 600

plt.rc("figure", dpi=dpi)

def fig_to_numpy(fig, **kwargs):
    io_buf = io.BytesIO()
    fig.savefig(io_buf, format="raw", **kwargs)
    io_buf.seek(0)
    img_arr = np.reshape(
        np.frombuffer(io_buf.getvalue(), dtype=np.uint8),
        newshape=(int(fig.bbox.bounds[3]), int(fig.bbox.bounds[2]), -1),
    )
    io_buf.close()
    return img_arr
colors = [(0.0, 'black'), (1.0, 'red')]

# Create a LinearSegmentedColormap object
def pure_colormap(color: str):
    return LinearSegmentedColormap.from_list(color, [(0, 0, 0, 1), color])

with plt.ioff():
    img = nib.load(snakemake.input["img"])
    vmax = np.mean(np.ma.masked_equal(np.sum(img.get_fdata(), axis=-1), 0)) * 2
    def get_layer(data, cmap):
        fig, axs = plt.subplots(3, 1, figsize=(30, 15), facecolor='black')
        for i in range(len(axs)):
            coords = np.r_[-50:50:7j]
            plotting.plot_img(
                data,
                coords,
                display_mode=["z", "x", "y"][i],
                black_bg=True,
                vmin=0,
                axes=axs[i],
                vmax=vmax,
                cmap=cmap,
            )
        plt.close(fig)
        w_ = 30 * dpi
        h_ = 15 * dpi
        crop_w = np.s_[int(h_ * 0.15) : int(h_ * 0.9), int(w_ * 0.15) : int(w_ * 0.9)]
        return fig_to_numpy(fig)[..., :3][crop_w]

    data = (
        get_layer(img.slicer[..., 0], pure_colormap('red')) + 
        get_layer(img.slicer[..., 1], pure_colormap('green')) + 
        get_layer(img.slicer[..., 2], pure_colormap('blue'))
    )
    fig, ax = plt.subplots(figsize=(30, 15))
    ax.imshow(data)
    ax.axes.axis('off')
    fig.savefig(snakemake.output.png, dpi=600, facecolor="black")
    plt.close(fig)
