
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from skimage.measure import label, regionprops_table
from skimage.filters import threshold_otsu
from skimage.morphology import remove_small_objects
import tifffile as tiff

def calculate_shape_descriptors_across_frames(stack, channel, min_size=5, organelle_name="organelle", output_dir="outputs"):
    all_shapes = []

    for f in range(stack.shape[0]):
        frame = stack[f, ..., channel]
        mask = frame > threshold_otsu(frame)
        mask = remove_small_objects(mask, min_size=min_size)
        labels = label(mask)

        props = regionprops_table(labels, properties=[
            "label", "area", "eccentricity", "perimeter", "solidity", "major_axis_length", "minor_axis_length"
        ])
        shape_df = pd.DataFrame(props)
        shape_df["Frame"] = f
        shape_df["Aspect_Ratio"] = shape_df["major_axis_length"] / shape_df["minor_axis_length"].replace(0, np.nan)
        shape_df["Circularity"] = (4 * np.pi * shape_df["area"]) / (shape_df["perimeter"]**2 + 1e-8)
        all_shapes.append(shape_df)

    shape_df_all = pd.concat(all_shapes, ignore_index=True)
    shape_df_all.to_csv(f"{output_dir}/{organelle_name}_shape_descriptors_all_frames.csv", index=False)

    plt.figure(figsize=(10, 5))
    sns.boxplot(data=shape_df_all, x="Frame", y="Aspect_Ratio")
    sns.lineplot(data=shape_df_all.groupby("Frame")["Aspect_Ratio"].mean().reset_index(), x="Frame", y="Aspect_Ratio", color="black", label="Mean Trend")
    plt.title(f"{organelle_name.capitalize()} Aspect Ratio Across Frames (with Trend)")
    plt.legend()
    plt.savefig(f"{output_dir}/{organelle_name}_aspect_ratio_trend.png")
    plt.close()

    plt.figure(figsize=(10, 5))
    sns.boxplot(data=shape_df_all, x="Frame", y="Circularity")
    sns.lineplot(data=shape_df_all.groupby("Frame")["Circularity"].mean().reset_index(), x="Frame", y="Circularity", color="black", label="Mean Trend")
    plt.title(f"{organelle_name.capitalize()} Circularity Across Frames (with Trend)")
    plt.legend()
    plt.savefig(f"{output_dir}/{organelle_name}_circularity_trend.png")
    plt.close()

# --- July 2025 Update ---
import pandas as pd
import matplotlib.pyplot as plt

def compute_and_export_shape_metrics(label_stack, save_prefix="shape_metrics"):
    from skimage.measure import regionprops
    all_data = []

    for t, frame in enumerate(label_stack):
        for region in regionprops(frame):
            label_id = region.label
            area = region.area
            perimeter = region.perimeter if region.perimeter > 0 else 1
            circularity = (4 * 3.14159 * area) / (perimeter ** 2)
            major = region.major_axis_length if region.major_axis_length > 0 else 1
            minor = region.minor_axis_length if region.minor_axis_length > 0 else 1
            aspect_ratio = major / minor
            all_data.append({
                "frame": t,
                "label": label_id,
                "area": area,
                "circularity": circularity,
                "aspect_ratio": aspect_ratio
            })

    df = pd.DataFrame(all_data)
    df.to_csv(f"{save_prefix}.csv", index=False)

    # Plot average circularity and aspect ratio over time
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    df.groupby("frame")["circularity"].mean().plot(ax=axes[0], title="Avg Circularity")
    df.groupby("frame")["aspect_ratio"].mean().plot(ax=axes[1], title="Avg Aspect Ratio")
    for ax in axes:
        ax.set_xlabel("Frame")
    plt.tight_layout()
    plt.savefig(f"{save_prefix}_plots.png")
    plt.close()

def run_shape_analysis_for_both_stacks(mito_stack, lyso_stack):
    compute_and_export_shape_metrics_with_trendline(mito_stack, "mitochondria_shape_metrics")
    compute_and_export_shape_metrics_with_trendline(lyso_stack, "lysosome_shape_metrics")
