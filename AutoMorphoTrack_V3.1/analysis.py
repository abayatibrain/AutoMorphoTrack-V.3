
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from skimage.measure import regionprops
from scipy.spatial.distance import cdist

def compute_shape_descriptors(data, output_dir):
    mito_labels = data['mito_labels']
    records = []

    for frame_idx, labels in enumerate(mito_labels):
        props = regionprops(labels)
        for p in props:
            area = p.area
            ecc = p.eccentricity
            perim = p.perimeter if p.perimeter > 0 else 1
            circ = (4 * np.pi * area) / (perim ** 2)
            solidity = p.solidity
            aspect_ratio = p.major_axis_length / p.minor_axis_length if p.minor_axis_length > 0 else 0
            records.append({
                "Frame": frame_idx,
                "Label": p.label,
                "Area": area,
                "Eccentricity": ecc,
                "Circularity": circ,
                "Solidity": solidity,
                "Aspect_Ratio": aspect_ratio
            })

    df = pd.DataFrame(records)
    df.to_csv(f"{output_dir}/shape_descriptors.csv", index=False)

    # Boxplots for visual analysis
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=df, x="Frame", y="Aspect_Ratio")
    plt.title("Mitochondrial Aspect Ratio per Frame")
    plt.savefig(f"{output_dir}/aspect_ratio_boxplot.png")
    plt.close()

    plt.figure(figsize=(10, 6))
    sns.boxplot(data=df, x="Frame", y="Circularity")
    plt.title("Mitochondrial Circularity per Frame")
    plt.savefig(f"{output_dir}/circularity_boxplot.png")
    plt.close()

def compute_directionality(tracking_csv_path, output_csv_path, output_img_path, frame_img, channel):
    df = pd.read_csv(tracking_csv_path)
    directionality = []
    grouped = df.groupby("Object_ID")

    for obj_id, group in grouped:
        group = group.sort_values("Frame")
        coords = group[["X", "Y"]].to_numpy()
        if len(coords) < 2:
            continue
        net_disp = np.linalg.norm(coords[-1] - coords[0])
        total_path = np.sum(np.linalg.norm(np.diff(coords, axis=0), axis=1))
        if total_path == 0:
            continue
        d_index = net_disp / total_path
        directionality.append({"Object_ID": obj_id, "Directionality": d_index})

    dir_df = pd.DataFrame(directionality)
    dir_df.to_csv(output_csv_path, index=False)

    # Scatter of directionality
    plt.figure(figsize=(8, 5))
    sns.histplot(dir_df["Directionality"], bins=20)
    plt.title("Distribution of Directionality Index")
    plt.savefig(output_img_path)
    plt.close()

def compute_lyso_mito_distances(data, output_dir):
    dist_records = []

    for frame_idx in range(len(data["mito_labels"])):
        mito_props = regionprops(data["mito_labels"][frame_idx])
        lyso_props = regionprops(data["lyso_labels"][frame_idx])
        if not mito_props or not lyso_props:
            continue

        mito_coords = np.array([p.centroid for p in mito_props])
        lyso_coords = np.array([p.centroid for p in lyso_props])
        dists = cdist(lyso_coords, mito_coords)
        min_dists = dists.min(axis=1)

        for d in min_dists:
            dist_records.append({"Frame": frame_idx, "MinDistance": d})

    df = pd.DataFrame(dist_records)
    df.to_csv(f"{output_dir}/lyso_mito_min_distances.csv", index=False)

    plt.figure(figsize=(8, 5))
    sns.histplot(df["MinDistance"], bins=30)
    plt.title("Minimum Lysosome-Mitochondrion Distances")
    plt.xlabel("Distance (pixels)")
    plt.savefig(f"{output_dir}/lyso_mito_distance_hist.png")
    plt.close()

# --- July 2025 Update: Colocalization Analysis ---
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from skimage.filters import threshold_otsu

def run_masked_pearson_colocalization(image_stack, save_prefix="pearson_colocalization"):
    pearson_results = []

    for t in range(image_stack.shape[0]):
        lyso = image_stack[t, 0].astype(float)
        mito = image_stack[t, 1].astype(float)

        mask_lyso = lyso > threshold_otsu(lyso)
        mask_mito = mito > threshold_otsu(mito)
        combined_mask = mask_lyso | mask_mito

        lyso_masked = lyso[combined_mask]
        mito_masked = mito[combined_mask]

        if len(lyso_masked) > 0 and len(mito_masked) > 0:
            r, _ = pearsonr(lyso_masked, mito_masked)
            pearson_results.append(abs(r))
        else:
            pearson_results.append(np.nan)

    df = pd.DataFrame({"frame": list(range(len(pearson_results))), "pearson_r": pearson_results})
    df.to_csv(f"{save_prefix}_per_frame.csv", index=False)

    plt.figure(figsize=(6, 4))
    plt.plot(df["frame"], df["pearson_r"], marker='o', color='c')
    plt.title("Pearson Colocalization (Masked) Across Frames")
    plt.xlabel("Frame")
    plt.ylabel("Pearson r")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"{save_prefix}_plot.png")
    plt.close()

    return df

def visualize_colocalization_overlay(image_stack, save_prefix="colocalization_overlay"):
    from skimage.filters import threshold_otsu

    selected_frames = np.linspace(0, image_stack.shape[0] - 1, 5, dtype=int)
    for t in selected_frames:
        lyso = image_stack[t, 0]
        mito = image_stack[t, 1]

        mask_lyso = lyso > threshold_otsu(lyso)
        mask_mito = mito > threshold_otsu(mito)
        coloc = mask_lyso & mask_mito

        overlay = np.zeros((*coloc.shape, 3), dtype=np.uint8)
        overlay[..., 0] = mito.astype(np.uint8)  # red
        overlay[..., 1] = lyso.astype(np.uint8)  # green
        overlay[coloc, :] = [0, 255, 255]  # cyan colocalization

        plt.imsave(f"{save_prefix}_frame{t}.png", overlay)

from PIL import Image
import os

def generate_colocalization_montage(selected_frame_paths, save_path="colocalization_overlay_montage.png"):
    images = [Image.open(path) for path in selected_frame_paths]
    widths, heights = zip(*(img.size for img in images))
    total_width = sum(widths)
    max_height = max(heights)
    montage = Image.new('RGB', (total_width, max_height))
    x_offset = 0
    for img in images:
        montage.paste(img, (x_offset, 0))
        x_offset += img.width
    montage.save(save_path)


import tifffile
import numpy as np
import pandas as pd
import os
from skimage.filters import threshold_isodata
from skimage.measure import label, regionprops
from skimage.morphology import remove_small_objects
from scipy.spatial import distance

def threshold_and_detect(image_stack, min_size=5):
    thresh_stack = np.array([frame > threshold_isodata(frame) for frame in image_stack])
    masks = [label(remove_small_objects(frame, min_size=min_size)) for frame in thresh_stack]
    props_list = [regionprops(mask) for mask in masks]
    return props_list, masks

def classify_mitochondria_morphology(props_list):
    morphology_per_frame = []
    for i, props in enumerate(props_list):
        elongated = 0
        punctate = 0
        for prop in props:
            area = prop.area
            eccentricity = prop.eccentricity
            if area >= 0.025 and eccentricity >= 0.85:
                elongated += 1
            elif area >= 0.03 and eccentricity <= 0.7:
                punctate += 1
        morphology_per_frame.append({'Frame': i, 'Elongated': elongated, 'Punctate': punctate})
    return pd.DataFrame(morphology_per_frame)

def calculate_motion(props_list):
    tracks = [[prop.centroid for prop in props] for props in props_list]
    displacements = []
    velocities = []
    for i in range(1, len(tracks)):
        prev_centroids = tracks[i - 1]
        curr_centroids = tracks[i]
        frame_displacements = []
        for curr in curr_centroids:
            if not prev_centroids:
                frame_displacements.append(0)
                continue
            dists = [distance.euclidean(curr, prev) for prev in prev_centroids]
            frame_displacements.append(np.min(dists))
        displacements.append(np.mean(frame_displacements) if frame_displacements else 0)
        velocities.append(np.mean(frame_displacements) if frame_displacements else 0)
    return pd.DataFrame({'Frame': range(1, len(displacements) + 1),
                         'Mean Displacement': displacements,
                         'Mean Velocity': velocities})

def run_full_analysis(image_path, output_dir, mito_channel=0, lyso_channel=1):
    os.makedirs(output_dir, exist_ok=True)

    img_stack = tifffile.imread(image_path)
    mito_stack = img_stack[:, mito_channel, :, :]
    lyso_stack = img_stack[:, lyso_channel, :, :]

    mito_props_list, _ = threshold_and_detect(mito_stack, min_size=30)
    lyso_props_list, _ = threshold_and_detect(lyso_stack, min_size=15)

    morphology_df = classify_mitochondria_morphology(mito_props_list)
    morphology_df.to_csv(os.path.join(output_dir, "mitochondria_morphology.csv"), index=False)

    mito_motion_df = calculate_motion(mito_props_list)
    lyso_motion_df = calculate_motion(lyso_props_list)

    mito_motion_df.to_csv(os.path.join(output_dir, "mitochondria_motion.csv"), index=False)
    lyso_motion_df.to_csv(os.path.join(output_dir, "lysosome_motion.csv"), index=False)


import tifffile
import numpy as np
import pandas as pd
import os
from scipy.spatial import distance
from skimage.measure import regionprops

from detection import detect_and_visualize

def classify_mitochondria_morphology(mito_labels):
    morphology_per_frame = []
    for i, labels in enumerate(mito_labels):
        props = regionprops(labels)
        elongated = 0
        punctate = 0
        for prop in props:
            area = prop.area
            eccentricity = prop.eccentricity
            if area >= 0.025 and eccentricity >= 0.85:
                elongated += 1
            elif area >= 0.03 and eccentricity <= 0.7:
                punctate += 1
        morphology_per_frame.append({'Frame': i, 'Elongated': elongated, 'Punctate': punctate})
    return pd.DataFrame(morphology_per_frame)

def calculate_motion(label_stack):
    displacements = []
    velocities = []
    for i in range(1, len(label_stack)):
        prev_props = regionprops(label_stack[i - 1])
        curr_props = regionprops(label_stack[i])
        prev_centroids = [p.centroid for p in prev_props]
        curr_centroids = [p.centroid for p in curr_props]
        frame_displacements = []
        for curr in curr_centroids:
            if not prev_centroids:
                frame_displacements.append(0)
                continue
            dists = [distance.euclidean(curr, prev) for prev in prev_centroids]
            frame_displacements.append(np.min(dists))
        displacements.append(np.mean(frame_displacements) if frame_displacements else 0)
        velocities.append(np.mean(frame_displacements) if frame_displacements else 0)
    return pd.DataFrame({'Frame': range(1, len(displacements) + 1),
                         'Mean Displacement': displacements,
                         'Mean Velocity': velocities})

def run_full_analysis(image_path, output_dir, mito_channel=0, lyso_channel=1):
    os.makedirs(output_dir, exist_ok=True)

    stack = tifffile.imread(image_path)
    n_frames = stack.shape[0]
    lyso_labels_stack = []
    mito_labels_stack = []

    for i in range(n_frames):
        frame = np.stack([stack[i, lyso_channel], stack[i, mito_channel]])
        lyso_labels, mito_labels = detect_and_visualize(frame, min_size=15, output_dir=output_dir)
        lyso_labels_stack.append(lyso_labels)
        mito_labels_stack.append(mito_labels)

    morphology_df = classify_mitochondria_morphology(mito_labels_stack)
    mito_motion_df = calculate_motion(mito_labels_stack)
    lyso_motion_df = calculate_motion(lyso_labels_stack)

    morphology_df.to_csv(os.path.join(output_dir, "mitochondria_morphology.csv"), index=False)
    mito_motion_df.to_csv(os.path.join(output_dir, "mitochondria_motion.csv"), index=False)
    lyso_motion_df.to_csv(os.path.join(output_dir, "lysosome_motion.csv"), index=False)


# --- New step: Lysosomal Count Analysis ---
from Lyso_count import count_lysosomes
from visualization import overlay_lysosome_count
import matplotlib.pyplot as plt

def run_lysosomal_count_analysis(lysosome_stack, lysosome_masks, output_dir):
    # Count lysosomes per frame and save CSV
    lyso_counts = count_lysosomes(lysosome_masks, output_dir)

    # Overlay counts on Frame 0
    overlay_lysosome_count(lysosome_stack[0], lysosome_masks[0], output_dir)

    # Plot count over time
    plt.figure()
    plt.plot(lyso_counts["Frame"], lyso_counts["Lysosome_Count"], marker='o')
    plt.xlabel("Frame")
    plt.ylabel("Lysosome Count")
    plt.title("Lysosome Count Over Time")
    plt.savefig(f"{output_dir}/lysosome_count_over_time.png", dpi=300)
    plt.close()



import matplotlib.pyplot as plt


import matplotlib.pyplot as plt


import matplotlib.pyplot as plt

def export_morphology_summary(df_morph, output_dir="outputs"):
    """
    Save mitochondrial morphology CSV + summary plots with percentages and stacked bar.
    Layout: 3x2 grid
      A: Counts (elongated vs punctate)
      B: Percentages (elongated vs punctate)
      C: Histogram of elongated counts
      D: Histogram of punctate counts
      E: Stacked bar plot of proportions per frame
      F: Empty
    """
    # Add percentage columns
    df_morph = df_morph.copy()
    df_morph["total"] = df_morph["elongated"] + df_morph["punctate"]
    df_morph["pct_elongated"] = df_morph.apply(lambda r: (r["elongated"]/r["total"]*100) if r["total"]>0 else 0, axis=1)
    df_morph["pct_punctate"] = df_morph.apply(lambda r: (r["punctate"]/r["total"]*100) if r["total"]>0 else 0, axis=1)

    # Save CSV with percentages
    csv_path = os.path.join(output_dir, "mitochondria_morphology.csv")
    df_morph.to_csv(csv_path, index=False)

    # Create figure with 3x2 grid (6 panels)
    fig, axes = plt.subplots(3, 2, figsize=(12, 15))

    # Panel A: Counts
    axes[0, 0].plot(df_morph["frame"], df_morph["elongated"], label="Elongated", color="blue")
    axes[0, 0].plot(df_morph["frame"], df_morph["punctate"], label="Punctate", color="green")
    axes[0, 0].set_xlabel("Frame"); axes[0, 0].set_ylabel("Count")
    axes[0, 0].set_title("A", loc="left", fontsize=14, fontweight="bold")
    axes[0, 0].legend()

    # Panel B: Percentages
    axes[0, 1].plot(df_morph["frame"], df_morph["pct_elongated"], label="% Elongated", color="blue")
    axes[0, 1].plot(df_morph["frame"], df_morph["pct_punctate"], label="% Punctate", color="green")
    axes[0, 1].set_xlabel("Frame"); axes[0, 1].set_ylabel("Percentage (%)")
    axes[0, 1].set_ylim(0, 100)
    axes[0, 1].set_title("B", loc="left", fontsize=14, fontweight="bold")
    axes[0, 1].legend()

    # Panel C: Histogram of elongated counts
    axes[1, 0].hist(df_morph["elongated"], bins=20, alpha=0.6, color="blue", label="Elongated")
    axes[1, 0].set_xlabel("Count"); axes[1, 0].set_ylabel("Frequency")
    axes[1, 0].set_title("C", loc="left", fontsize=14, fontweight="bold")
    axes[1, 0].legend()

    # Panel D: Histogram of punctate counts
    axes[1, 1].hist(df_morph["punctate"], bins=20, alpha=0.6, color="green", label="Punctate")
    axes[1, 1].set_xlabel("Count"); axes[1, 1].set_ylabel("Frequency")
    axes[1, 1].set_title("D", loc="left", fontsize=14, fontweight="bold")
    axes[1, 1].legend()

    # Panel E: Stacked bar plot (elongated vs punctate proportions)
    axes[2, 0].bar(df_morph["frame"], df_morph["elongated"], color="blue", label="Elongated")
    axes[2, 0].bar(df_morph["frame"], df_morph["punctate"], bottom=df_morph["elongated"], color="green", label="Punctate")
    axes[2, 0].set_xlabel("Frame"); axes[2, 0].set_ylabel("Count")
    axes[2, 0].set_title("E", loc="left", fontsize=14, fontweight="bold")
    axes[2, 0].legend()

    # Panel F: Empty
    axes[2, 1].axis("off")
    axes[2, 1].set_title("F", loc="left", fontsize=14, fontweight="bold")

    plt.tight_layout()
    fig_path = os.path.join(output_dir, "morphology_summary.png")
    plt.savefig(fig_path, dpi=150)
    plt.close(fig)
    print(f"Saved morphology summary: {fig_path}")
def export_lysosome_summary(df_lyso, output_dir="outputs"):
    """
    Save lysosome counts CSV + summary plots.
    - Line plot: lysosome counts per frame
    - Histogram: lysosome counts distribution
    """
    csv_path = os.path.join(output_dir, "lysosome_counts.csv")
    df_lyso.to_csv(csv_path, index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Line plot
    axes[0].plot(df_lyso["frame"], df_lyso["count"], color="green")
    axes[0].set_xlabel("Frame"); axes[0].set_ylabel("Count")
    axes[0].set_title("Lysosome Counts per Frame")

    # Histogram
    axes[1].hist(df_lyso["count"], bins=20, alpha=0.7, color="green")
    axes[1].set_xlabel("Count"); axes[1].set_ylabel("Frequency")
    axes[1].set_title("Distribution of Lysosome Counts")

    plt.tight_layout()
    fig_path = os.path.join(output_dir, "lysosome_summary.png")
    plt.savefig(fig_path, dpi=150)
    plt.close(fig)
    print(f"Saved lysosome summary: {fig_path}")



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import cv2
from skimage.measure import label, regionprops

def compute_masked_pearson(mito_img, lyso_img):
    """
    Compute masked Pearson correlation between mitochondria and lysosomes.
    Mask = union of thresholded mito + lyso.
    """
    mito_norm = cv2.normalize(mito_img, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    lyso_norm = cv2.normalize(lyso_img, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)

    _, mito_mask = cv2.threshold(mito_norm, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    _, lyso_mask = cv2.threshold(lyso_norm, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)

    mask = ((mito_mask > 0) | (lyso_mask > 0))
    if np.sum(mask) == 0:
        return 0

    mito_vals = mito_norm[mask].astype(np.float32)
    lyso_vals = lyso_norm[mask].astype(np.float32)

    if len(mito_vals) < 2 or len(lyso_vals) < 2:
        return 0

    r, _ = pearsonr(mito_vals, lyso_vals)
    return r

def export_colocalization(stack, output_dir="outputs"):
    """
    Compute masked Pearson correlation for all frames and export results.
    - Saves CSV (colocalization.csv)
    - Saves summary figure (line plot + histogram)
    """
    coloc_results = []
    for f in range(stack.shape[0]):
        r = compute_masked_pearson(stack[f,0], stack[f,1])
        coloc_results.append({"frame": f, "pearson_r": r})

    df_coloc = pd.DataFrame(coloc_results)
    csv_path = os.path.join(output_dir, "colocalization.csv")
    df_coloc.to_csv(csv_path, index=False)

    # Create summary figure
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    axes[0].plot(df_coloc["frame"], df_coloc["pearson_r"], color="purple")
    axes[0].set_xlabel("Frame"); axes[0].set_ylabel("Pearson r")
    axes[0].set_title("A", loc="left", fontsize=14, fontweight="bold")

    axes[1].hist(df_coloc["pearson_r"], bins=20, color="purple", alpha=0.7)
    axes[1].set_xlabel("Pearson r"); axes[1].set_ylabel("Frequency")
    axes[1].set_title("B", loc="left", fontsize=14, fontweight="bold")

    plt.tight_layout()
    fig_path = os.path.join(output_dir, "colocalization_summary.png")
    plt.savefig(fig_path, dpi=150)
    plt.close(fig)
    print(f"Saved colocalization summary: {fig_path}")

def visualize_colocalization_frame(stack, frame_idx=0):
    """
    Create overlay visualization of colocalization for a single frame.
    Mito = red, Lyso = green, Colocalized = yellow
    """
    mito = stack[frame_idx,0]
    lyso = stack[frame_idx,1]
    mito_norm = cv2.normalize(mito, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    lyso_norm = cv2.normalize(lyso, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)

    _, mito_mask = cv2.threshold(mito_norm, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    _, lyso_mask = cv2.threshold(lyso_norm, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)

    overlay = np.zeros((mito.shape[0], mito.shape[1], 3), dtype=np.uint8)
    overlay[...,0] = mito_norm
    overlay[...,1] = lyso_norm

    coloc = (mito_mask > 0) & (lyso_mask > 0)
    overlay[coloc] = [255, 255, 0]  # yellow for overlap

    plt.figure(figsize=(6,6))
    plt.imshow(overlay)
    plt.axis("off")
    plt.show()



import os
import cv2

def images_to_video(image_dir, output_path, fps=5):
    """
    Convert sequence of saved .png images to MP4 video.
    """
    import glob
    images = sorted(glob.glob(os.path.join(image_dir, "*.png")))
    if not images:
        return
    frame = cv2.imread(images[0])
    height, width, _ = frame.shape
    fourcc = cv2.VideoWriter_fourcc(*"mp4v")
    video = cv2.VideoWriter(output_path, fourcc, fps, (width, height))
    for img_path in images:
        frame = cv2.imread(img_path)
        video.write(frame)
    video.release()

def run_full_pipeline(stack, output_dir="outputs"):
    """
    Run the full AutoMorphoTrack pipeline with PNG + MP4 outputs at each step.
    Steps: Thresholding -> Detection -> Morphology -> Lysosome counts -> Displacement -> Colocalization
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Thresholding video
    thr_dir = os.path.join(output_dir, "threshold")
    if os.path.exists(thr_dir):
        images_to_video(thr_dir, os.path.join(output_dir, "threshold_overlay.mp4"), fps=5)

    # Detection video
    det_dir = os.path.join(output_dir, "detection")
    if os.path.exists(det_dir):
        images_to_video(det_dir, os.path.join(output_dir, "detection_overlay.mp4"), fps=5)

    # Morphology video
    morph_dir = os.path.join(output_dir, "morphology")
    if os.path.exists(morph_dir):
        images_to_video(morph_dir, os.path.join(output_dir, "morphology_overlay.mp4"), fps=5)

    # Lysosome counts video
    lyso_dir = os.path.join(output_dir, "lysosomes")
    if os.path.exists(lyso_dir):
        images_to_video(lyso_dir, os.path.join(output_dir, "lysosomes.mp4"), fps=5)

    # Displacement videos
    disp_dir = os.path.join(output_dir, "displacement")
    if os.path.exists(disp_dir):
        images_to_video(os.path.join(disp_dir, "mitochondria"),
                        os.path.join(output_dir, "displacement_mitochondria.mp4"), fps=5)
        images_to_video(os.path.join(disp_dir, "lysosomes"),
                        os.path.join(output_dir, "displacement_lysosomes.mp4"), fps=5)
        images_to_video(os.path.join(disp_dir, "combined"),
                        os.path.join(output_dir, "displacement_combined.mp4"), fps=5)

    # Colocalization video
    coloc_dir = os.path.join(output_dir, "colocalization")
    if os.path.exists(coloc_dir):
        images_to_video(coloc_dir, os.path.join(output_dir, "colocalization_overlay.mp4"), fps=5)

    print(f"Full pipeline completed. Outputs saved to {output_dir}")
