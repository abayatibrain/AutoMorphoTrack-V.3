
import numpy as np
import matplotlib.pyplot as plt
from skimage.filters import threshold_otsu
from skimage.morphology import remove_small_objects
from skimage.measure import label
from skimage.segmentation import find_boundaries

def detect_and_visualize(frame, min_size=5, output_dir="outputs"):
    lyso = frame[0, :, :]
    mito = frame[1, :, :]

    lyso_mask = lyso > threshold_otsu(lyso)
    lyso_mask = remove_small_objects(lyso_mask, min_size=min_size)
    lyso_labels = label(lyso_mask)
    lyso_outline = find_boundaries(lyso_mask, mode='outer')

    mito_mask = mito > threshold_otsu(mito)
    mito_mask = remove_small_objects(mito_mask, min_size=min_size)
    mito_labels = label(mito_mask)
    mito_outline = find_boundaries(mito_mask, mode='outer')

    # Grayscale overlays with red outlines
    lyso_rgb = np.stack([lyso]*3, axis=-1)
    mito_rgb = np.stack([mito]*3, axis=-1)
    lyso_rgb[lyso_outline, 0] = 255
    lyso_rgb[lyso_outline, 1] = 0
    lyso_rgb[lyso_outline, 2] = 0
    mito_rgb[mito_outline, 0] = 255
    mito_rgb[mito_outline, 1] = 0
    mito_rgb[mito_outline, 2] = 0

    # Save side-by-side
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    axes[0].imshow(lyso_rgb.astype(np.uint8))
    axes[0].set_title("Lysosomes with Red Outlines")
    axes[0].axis("off")

    axes[1].imshow(mito_rgb.astype(np.uint8))
    axes[1].set_title("Mitochondria with Red Outlines")
    axes[1].axis("off")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/individual_channels_with_outlines.png")
    plt.close()

    return lyso_labels, mito_labels



def detect_and_visualize(frame, min_size=5, output_dir="outputs"):
    import numpy as np
    import matplotlib.pyplot as plt
    from skimage.filters import threshold_otsu
    from skimage.morphology import remove_small_objects
    from skimage.measure import label
    from skimage.segmentation import find_boundaries
    import os
    from PIL import Image

    lyso = frame[0, :, :]
    mito = frame[1, :, :]

    # Save individual grayscale channels
    Image.fromarray(lyso.astype(np.uint8)).save(os.path.join(output_dir, "Channel_0_Lysosomes_Frame0.png"))
    Image.fromarray(mito.astype(np.uint8)).save(os.path.join(output_dir, "Channel_1_Mitochondria_Frame0.png"))

    # Process and outline
    lyso_mask = lyso > threshold_otsu(lyso)
    lyso_mask = remove_small_objects(lyso_mask, min_size=min_size)
    lyso_labels = label(lyso_mask)
    lyso_outline = find_boundaries(lyso_mask, mode='outer')

    mito_mask = mito > threshold_otsu(mito)
    mito_mask = remove_small_objects(mito_mask, min_size=min_size)
    mito_labels = label(mito_mask)
    mito_outline = find_boundaries(mito_mask, mode='outer')

    lyso_rgb = np.stack([lyso]*3, axis=-1)
    mito_rgb = np.stack([mito]*3, axis=-1)
    lyso_rgb[lyso_outline, 0] = 255
    lyso_rgb[lyso_outline, 1] = 0
    lyso_rgb[lyso_outline, 2] = 0
    mito_rgb[mito_outline, 0] = 255
    mito_rgb[mito_outline, 1] = 0
    mito_rgb[mito_outline, 2] = 0

    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    axes[0].imshow(lyso_rgb.astype(np.uint8))
    axes[0].axis("off")
    axes[1].imshow(mito_rgb.astype(np.uint8))
    axes[1].axis("off")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "Detection_Outlines_Frame0.png"))
    plt.close()

    return lyso_labels, mito_labels



def count_lysosomes_across_frames(image_stack, output_dir, min_size=10):
    import os
    import pandas as pd
    import matplotlib.pyplot as plt
    from skimage.filters import threshold_otsu
    from skimage.morphology import remove_small_objects
    from skimage.measure import label, regionprops
    from skimage.segmentation import find_boundaries
    import numpy as np

    os.makedirs(output_dir, exist_ok=True)
    counts = []

    for i, image in enumerate(image_stack):
        thresh = threshold_otsu(image)
        mask = image > thresh
        mask = remove_small_objects(mask, min_size=min_size)
        labels = label(mask)
        count = len(regionprops(labels))
        counts.append({"Frame": i, "Lysosome_Count": count})

        # Save labeled overlay image
        outline = find_boundaries(mask, mode='outer')
        overlay = np.stack([image]*3, axis=-1)
        overlay[outline, 1] = 255  # Green outlines
        overlay[outline, 0] = 0
        overlay[outline, 2] = 0

        plt.figure(figsize=(5, 5))
        plt.imshow(overlay.astype(np.uint8))
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"Lysosome_Count_Frame{i}.png"))
        plt.close()

    df = pd.DataFrame(counts)
    df.to_csv(os.path.join(output_dir, "lysosome_counts.csv"), index=False)

    # Plot count
    plt.figure(figsize=(8, 4))
    plt.plot(df["Frame"], df["Lysosome_Count"], marker="o", color="green")
    plt.xlabel("Frame")
    plt.ylabel("Lysosome Count")
    plt.title("Lysosome Count Across Frames")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "Lysosome_Counts_Plot.png"))
    plt.close()
