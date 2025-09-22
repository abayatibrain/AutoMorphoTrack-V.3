
import numpy as np
import pandas as pd
from skimage.measure import regionprops
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap

def track_organelles(data, output_dir):
    def extract_tracks(labels):
        tracks = {}
        for frame_idx, frame in enumerate(labels):
            props = regionprops(frame)
            for obj in props:
                label_id = obj.label
                if label_id not in tracks:
                    tracks[label_id] = []
                tracks[label_id].append((frame_idx, obj.centroid[0], obj.centroid[1]))
        return tracks

    mito_tracks = extract_tracks(data['mito_labels'])
    lyso_tracks = extract_tracks(data['lyso_labels'])

    # Flatten for CSV export
    mito_data = [{'Object_ID': k, 'Frame': f, 'Y': y, 'X': x}
                 for k, v in mito_tracks.items() for f, y, x in v]
    lyso_data = [{'Object_ID': k, 'Frame': f, 'Y': y, 'X': x}
                 for k, v in lyso_tracks.items() for f, y, x in v]

    pd.DataFrame(mito_data).to_csv(f"{output_dir}/displacement_velocity_mito.csv", index=False)
    pd.DataFrame(lyso_data).to_csv(f"{output_dir}/displacement_velocity_lyso.csv", index=False)

    visualize_multicolor_tracks(data, mito_tracks, lyso_tracks, output_dir)

def visualize_multicolor_tracks(data, mito_tracks, lyso_tracks, output_dir):
    cmap = get_cmap("tab20")
    raw = data['stack'][0]

    def draw_tracks(ax, tracks, channel, color_shift=0):
        ax.imshow(raw[..., channel], cmap='gray')
        for i, (obj_id, points) in enumerate(tracks.items()):
            if len(points) < 2:
                continue
            coords = np.array([(x, y) for _, y, x in points])
            ax.plot(coords[:, 0], coords[:, 1], color=cmap((i + color_shift) % 20), linewidth=1)

        ax.axis('off')

    # Lysosomal tracks
    fig1, ax1 = plt.subplots()
    draw_tracks(ax1, lyso_tracks, channel=0)
    fig1.savefig(f"{output_dir}/tracklines_lyso_frame0.png", bbox_inches='tight', pad_inches=0)
    plt.close(fig1)

    # Mitochondrial tracks
    fig2, ax2 = plt.subplots()
    draw_tracks(ax2, mito_tracks, channel=1, color_shift=10)
    fig2.savefig(f"{output_dir}/tracklines_mito_frame0.png", bbox_inches='tight', pad_inches=0)
    plt.close(fig2)


# --- July 2025 Update ---
from skimage.measure import label as sk_label
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def detect_organelles(image, min_size=30):
    from skimage.filters import threshold_otsu
    from skimage.morphology import remove_small_objects, binary_opening, disk
    thresh = threshold_otsu(image)
    binary = image > thresh
    cleaned = binary_opening(binary, disk(1))
    cleaned = remove_small_objects(cleaned, min_size=min_size)
    labeled = sk_label(cleaned)
    return labeled

def track_centroids(label_stack):
    from skimage.measure import regionprops
    tracks = {}
    for t, frame in enumerate(label_stack):
        for region in regionprops(frame):
            label_id = region.label
            y, x = region.centroid
            if label_id not in tracks:
                tracks[label_id] = []
            tracks[label_id].append((t, x, y))
    return tracks

def plot_tracks(frame, tracks, title, save_path=None):
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.imshow(frame, cmap='gray')
    colormap = cm.get_cmap('nipy_spectral', len(tracks))
    for i, (label, points) in enumerate(tracks.items()):
        color = colormap(i)
        for t, x, y in points:
            ax.plot(x, y, 'o', markersize=2, color=color)
    ax.set_title(title)
    ax.axis('off')
    if save_path:
        plt.savefig(save_path, bbox_inches='tight')
    plt.close()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import Normalize
from PIL import Image

def draw_arrowed_trajectory(frame, tracks, title, save_path=None):
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.imshow(frame, cmap='gray')
    cmap = get_cmap('nipy_spectral')
    for label_id, points in tracks.items():
        if len(points) < 2:
            continue
        x_coords = [x for _, x, _ in points]
        y_coords = [y for _, _, y in points]
        colors = cmap(np.linspace(0, 1, len(points)))
        for i in range(len(points) - 1):
            ax.annotate("",
                        xy=(x_coords[i+1], y_coords[i+1]),
                        xytext=(x_coords[i], y_coords[i]),
                        arrowprops=dict(arrowstyle="->", color=colors[i], lw=1),
                        annotation_clip=False)
    ax.set_title(title)
    ax.axis('off')
    if save_path:
        plt.savefig(save_path, bbox_inches='tight')
    plt.close()

def draw_rainbow_trails(frame, tracks, title, save_path=None):
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.imshow(frame, cmap='gray')
    cmap = get_cmap('plasma')
    for label_id, points in tracks.items():
        if len(points) < 2:
            continue
        x_coords = [x for _, x, _ in points]
        y_coords = [y for _, _, y in points]
        trail_colors = cmap(np.linspace(0, 1, len(points)))
        for i in range(len(points) - 1):
            ax.plot([x_coords[i], x_coords[i+1]],
                    [y_coords[i], y_coords[i+1]],
                    color=trail_colors[i], lw=1)
    ax.set_title(title)
    ax.axis('off')
    if save_path:
        plt.savefig(save_path, bbox_inches='tight')
    plt.close()

def draw_occupancy_heatmap(label_stack, title, save_path=None):
    occupancy = np.zeros(label_stack[0].shape, dtype=int)
    for labeled in label_stack:
        occupancy += labeled > 0
    plt.figure(figsize=(5, 5))
    plt.imshow(occupancy, cmap='hot')
    plt.title(title)
    plt.axis('off')
    if save_path:
        plt.savefig(save_path, bbox_inches='tight')
    plt.close()

def draw_start_end_overlay(frame, tracks, title, save_path=None):
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.imshow(frame, cmap='gray')
    for label_id, points in tracks.items():
        if len(points) < 2:
            continue
        _, x_start, y_start = points[0]
        _, x_end, y_end = points[-1]
        ax.plot([x_start, x_end], [y_start, y_end], color='yellow', lw=1)
        ax.plot(x_start, y_start, 'go', markersize=3)
        ax.plot(x_end, y_end, 'ro', markersize=3)
    ax.set_title(title)
    ax.axis('off')
    if save_path:
        plt.savefig(save_path, bbox_inches='tight')
    plt.close()
