import matplotlib.pyplot as plt

def plot_organelles(image, mask):
    plt.imshow(image, cmap='gray')
    plt.contour(mask, colors='red')
    plt.show()


import matplotlib.pyplot as plt
import numpy as np

def visualize_colocalization_frame0(mito_stack, lyso_stack):
    frame0_mito = mito_stack[0]
    frame0_lyso = lyso_stack[0]

    # Normalize images for display
    norm_mito = frame0_mito / np.max(frame0_mito)
    norm_lyso = frame0_lyso / np.max(frame0_lyso)

    # Colocalized pixels where both mito and lyso intensities are non-zero
    colocalized = np.logical_and(norm_mito > 0.1, norm_lyso > 0.1)

    rgb = np.zeros((*frame0_mito.shape, 3))
    rgb[..., 1] = norm_mito  # Red channel
    rgb[..., 0] = norm_lyso  # Green channel
    rgb[..., 2] = colocalized.astype(float)  # Cyan = Blue + Green, but we'll mark it with Blue

    plt.figure(figsize=(6, 6))
    plt.imshow(rgb)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig("frame0_colocalization_cyan.png", dpi=300, bbox_inches='tight')
    plt.close()


import matplotlib.pyplot as plt
import numpy as np

def draw_mitochondrial_morphology(frame_mask, props, filename="mitochondria_morphology_only.png"):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.imshow(frame_mask, cmap='Reds')
    for prop in props:
        y, x = prop.centroid
        if prop.area >= 0.2 and prop.eccentricity >= 0.85:
            ax.text(x, y, 'E', color='red', fontsize=12, ha='center', va='center')
        elif prop.area >= 0.03 and prop.eccentricity <= 0.7:
            ax.text(x, y, 'P', color='blue', fontsize=12, ha='center', va='center')
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

def draw_lysosome_count(frame_mask, props, filename="lysosome_count_only.png"):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.imshow(frame_mask, cmap='Greens')
    for i, prop in enumerate(props, 1):
        y, x = prop.centroid
        ax.text(x, y, str(i), color='lime', fontsize=8, ha='center', va='center')
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()


import matplotlib.pyplot as plt
from skimage.measure import label, regionprops

def overlay_lysosome_count(frame_img, mask, output_dir):
    """
    Overlay lysosome count numbers on Frame 0.
    """
    labeled = label(mask)
    props = regionprops(labeled)

    fig, ax = plt.subplots()
    ax.imshow(frame_img, cmap='gray')
    for i, prop in enumerate(props, start=1):
        y, x = prop.centroid
        ax.text(x, y, str(i), color='lime', fontsize=10, ha='center', va='center')

    ax.set_axis_off()
    save_path = f"{output_dir}/lysosome_count_overlay_Frame0.png"
    plt.savefig(save_path, bbox_inches='tight', dpi=300)
    plt.close(fig)



import matplotlib.pyplot as plt

def show_three_panel(overlay, panel_b, panel_c, labels=("A","B","C")):
    """
    Show three panels with overlay first, then individual channels.
    No titles, only A,B,C labels.
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    axes[0].imshow(overlay)
    axes[0].set_title(labels[0], loc="left", fontsize=14, fontweight="bold")
    axes[0].axis("off")

    axes[1].imshow(panel_b)
    axes[1].set_title(labels[1], loc="left", fontsize=14, fontweight="bold")
    axes[1].axis("off")

    axes[2].imshow(panel_c)
    axes[2].set_title(labels[2], loc="left", fontsize=14, fontweight="bold")
    axes[2].axis("off")

    plt.tight_layout()
    plt.show()
