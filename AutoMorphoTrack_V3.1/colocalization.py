
import os
import numpy as np
import matplotlib.pyplot as plt

def _normalize(img):
    img = img.astype(np.float32)
    vmin, vmax = np.nanmin(img), np.nanmax(img)
    if vmax <= vmin:
        return np.zeros_like(img, dtype=np.float32)
    return (img - vmin) / (vmax - vmin)

def masked_pearson(x, y, mask=None):
    """
    Compute Pearson correlation on pixels where mask==True.
    If mask is None, uses non-zero pixels in both x and y.
    Returns np.nan if < 3 valid pixels.
    """
    x = np.asarray(x)
    y = np.asarray(y)
    if mask is None:
        mask = (x > 0) & (y > 0)
    else:
        mask = mask.astype(bool)
    xv = x[mask].astype(np.float32)
    yv = y[mask].astype(np.float32)
    if xv.size < 3:
        return np.nan
    xv = xv - xv.mean()
    yv = yv - yv.mean()
    denom = (np.sqrt((xv**2).sum()) * np.sqrt((yv**2).sum()))
    if denom == 0:
        return np.nan
    return float((xv @ yv) / denom)

def compute_masked_pearson_per_frame(mito_stack, lyso_stack, mito_masks=None, lyso_masks=None):
    """
    mito_stack, lyso_stack: arrays of shape (T, H, W)
    Optionally provide mito_masks and lyso_masks (boolean stacks). If provided,
    the combined mask per frame is (mito_masks[t] | lyso_masks[t]) & (mito>0) & (lyso>0)
    Returns: list of r values per frame
    """
    mito_stack = np.asarray(mito_stack)
    lyso_stack = np.asarray(lyso_stack)
    T = min(mito_stack.shape[0], lyso_stack.shape[0])
    r_values = []
    for t in range(T):
        m = None
        if mito_masks is not None or lyso_masks is not None:
            mm = mito_masks[t] if mito_masks is not None else None
            lm = lyso_masks[t] if lyso_masks is not None else None
            if mm is not None and lm is not None:
                m = (mm | lm)
            elif mm is not None:
                m = mm
            elif lm is not None:
                m = lm
        r = masked_pearson(lyso_stack[t], mito_stack[t], mask=m)
        r_values.append(r)
    return r_values

def visualize_colocalization_map(frame_mito, frame_lyso, out_path="colocalization_frame0.png"):
    """
    Creates a side-by-side visualization:
    Left: RGB composite (R=mito, G=lyso) for the given frame.
    Right: Colocalization heatmap defined as normalized intensity product.
    """
    mito_n = _normalize(frame_mito)
    lyso_n = _normalize(frame_lyso)
    overlap = mito_n * lyso_n

    fig = plt.figure(figsize=(10, 4))
    # Left: composite
    ax1 = fig.add_subplot(1, 2, 1)
    composite = np.stack([mito_n, lyso_n, np.zeros_like(mito_n)], axis=-1)
    ax1.imshow(composite)
    ax1.set_title("Composite (R=mito, G=lyso)")
    ax1.axis("off")

    # Right: heatmap of overlap
    ax2 = fig.add_subplot(1, 2, 2)
    im = ax2.imshow(overlap, cmap="inferno")
    ax2.set_title("Colocalization heatmap (intensity product)")
    ax2.axis("off")
    cbar = plt.colorbar(im, ax=ax2, fraction=0.046, pad=0.04)
    cbar.set_label("Overlap intensity (arb. units)")

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path

def plot_correlation_timeseries(r_values, out_path="colocalization_timeseries.png"):
    """
    Line plot of masked Pearson r per frame.
    """
    x = np.arange(len(r_values))
    fig = plt.figure(figsize=(6, 4))
    plt.plot(x, r_values, marker="o")
    plt.xlabel("Frame")
    plt.ylabel("Masked Pearson r")
    plt.title("Mitochondria–Lysosome Colocalization Over Time")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path

def save_correlation_csv(r_values, out_csv="colocalization_timeseries.csv"):
    import pandas as pd
    df = pd.DataFrame({"frame": np.arange(len(r_values)), "masked_pearson_r": r_values})
    df.to_csv(out_csv, index=False)
    return out_csv

def run_colocalization_pipeline(mito_stack, lyso_stack, output_dir, mito_masks=None, lyso_masks=None):
    """
    Convenience function that:
    - computes masked Pearson per frame
    - saves CSV
    - saves line plot
    - saves a frame0 composite + heatmap figure
    Returns dictionary of output paths.
    """
    os.makedirs(output_dir, exist_ok=True)
    r_values = compute_masked_pearson_per_frame(mito_stack, lyso_stack, mito_masks, lyso_masks)
    csv_path = os.path.join(output_dir, "colocalization_timeseries.csv")
    png_ts = os.path.join(output_dir, "colocalization_timeseries.png")
    png_frame0 = os.path.join(output_dir, "colocalization_frame0.png")

    save_correlation_csv(r_values, csv_path)
    plot_correlation_timeseries(r_values, png_ts)
    # Use min length frame 0 safely
    f0_mito = mito_stack[0]
    f0_lyso = lyso_stack[0]
    visualize_colocalization_map(f0_mito, f0_lyso, png_frame0)

    return {"csv": csv_path, "timeseries_plot": png_ts, "frame0_coloc_map": png_frame0, "r_values": r_values}
