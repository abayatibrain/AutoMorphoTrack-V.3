
# Using AutoMorphoTrack on a Local Machine

To analyze organelle morphology, dynamics, and colocalization from time-lapse microscopy data, we developed and used the AutoMorphoTrack Python package. Below is a step-by-step guide on how to install and run the package locally, with details on expected outputs at each stage.

---

## 🔧 1. Installation

AutoMorphoTrack is compatible with Python ≥3.8. To install locally:

```bash
git clone https://github.com/yourusername/AutoMorphoTrack.git
cd AutoMorphoTrack
pip install -e .
```

Required dependencies:

```bash
pip install numpy scipy opencv-python matplotlib scikit-image pandas tifffile
```

---

## 🧬 2. Input Data Requirements

- Input should be a **2-channel time-lapse `.tif` stack**:
  - **Channel 0**: Lysosomes (punctate)
  - **Channel 1**: Mitochondria (filamentous)
- File format: `.tif` with shape `(frames, height, width, channels)`

---

## 🚀 3. Running the Analysis

Use the provided Jupyter notebook `AutoMorphoTrack_Full_Package.ipynb`, or run the pipeline via script:

```python
from AutoMorphoTrack import core

core.run_pipeline(
    filepath="path/to/your_image_stack.tif",
    output_dir="path/to/output_folder",
    lyso_channel=0,
    mito_channel=1,
    frame_to_visualize=0
)
```

This executes the full pipeline across all frames and generates visualizations on **Frame 0**.

---

## 📤 4. Outputs

The pipeline produces the following outputs in the specified `output_dir`:

### a. Detection and Thresholding

- `threshold_lyso_frame0.png` — Binary lysosome mask
- `threshold_mito_frame0.png` — Binary mitochondria mask
- `detection_frame0.png` — Organelle overlay (lysosomes = green, mitochondria = red)

### b. Morphology

- `morphology_frame0.png` — Each mitochondrion labeled as **Elongated (E)** or **Punctate (P)**
- `morphology_summary.csv` — Frame-by-frame classification table (area, eccentricity, class)

### c. Lysosome Quantification

- `lysosome_labels_frame0.png` — Individual lysosome IDs overlaid
- `lysosome_count.csv` — Count of lysosomes per frame

### d. Tracking

- Dot plots:
  - `tracks_lyso_frame0.png` — Lysosome positions across time
  - `tracks_mito_frame0.png` — Mitochondria positions across time
- Line trajectories:
  - `tracklines_lyso_frame0.png` — Lysosome tracks across frames, colored by object
  - `tracklines_mito_frame0.png` — Mitochondria tracks across frames, colored by object
- Tables:
  - `displacement_velocity_lyso.csv`
  - `displacement_velocity_mito.csv`

### e. Colocalization (if enabled)

- `colocalization_frame0.png` — Overlapping lysosome/mitochondria regions in **cyan**
- `pearson_colocalization.csv` — Framewise Pearson r values

---

## 🔄 5. Analysis Workflow

You can use each output to answer specific biological questions:

| Biological Question                  | Output                                 |
|--------------------------------------|----------------------------------------|
| How many lysosomes per neuron?       | `lysosome_count.csv`                   |
| Are mitochondria fragmented or tubular? | `morphology_summary.csv`, `morphology_frame0.png` |
| Do mitochondria move more under condition X? | `displacement_velocity_mito.csv` |
| Are mitochondria and lysosomes colocalized? | `colocalization_frame0.png`, `pearson_colocalization.csv` |
| Are dynamics altered over time?      | `tracklines_*.png`, trajectory CSVs    |

---

## 📈 6. Post-Processing and Quantification

- Visual outputs (e.g., `.png`) can be used directly for figures
- CSVs can be imported into Excel, GraphPad Prism, or Python/R for statistical comparison

---

## 📎 7. Reproducibility and Customization

- All detection thresholds and classification criteria are exposed in `core.py`
- Users can modify:
  - Morphology thresholds (area, eccentricity)
  - Colocalization metric (e.g., Manders, masked Pearson)
  - Track visualization style


---

## 📊 Additional Analyses

### 1. Mitochondrial Shape Descriptors
- Computes circularity, solidity, and aspect ratio
- Visualized with boxplots: `aspect_ratio_boxplot.png`, `circularity_boxplot.png`

### 2. Directionality Index
- Measures straightness of trajectories (net displacement / total path)
- Output: `directionality_summary.csv`, `directionality_hist.png`

### 3. Lysosome-Mitochondrion Proximity
- Calculates minimum distance of each lysosome to nearest mitochondrion per frame
- Output: `lyso_mito_min_distances.csv`, `lyso_mito_distance_hist.png`


## 🔄 Default Detection Criteria
AutoMorphoTrack uses relaxed detection thresholds by default:

- **Otsu thresholding** is applied to both channels
- **Small objects** (≥ 5 pixels) are retained for analysis
- This improves sensitivity for:
  - Small lysosomes
  - Filamentous mitochondria in dense clusters


### Shape Descriptor Analysis
AutoMorphoTrack includes a shape analysis module:

- `calculate_shape_descriptors_across_frames(stack, channel, ...)`
- Outputs per-organelle descriptors per frame
- CSVs and trendline plots for:
  - Aspect Ratio
  - Circularity


### Detection Visualization Update
AutoMorphoTrack now includes:
- Overlay of detected masks onto grayscale images
- Red outlines to help visualize organelle boundaries


## Colocalization: Analysis + Visualization

- **Function:** `run_colocalization_pipeline(mito_stack, lyso_stack, output_dir, mito_masks=None, lyso_masks=None)`
- **Outputs:**
  - `colocalization_timeseries.csv` — masked Pearson r per frame
  - `colocalization_timeseries.png` — line plot of r over time
  - `colocalization_frame0.png` — composite (R=mito, G=lyso) + heatmap of intensity product
- **Masks (optional):** Provide boolean stacks to restrict analysis to organelle-positive pixels.

**Example**
```python
from colocalization import run_colocalization_pipeline
out = run_colocalization_pipeline(mito_stack, lyso_stack, output_dir="outputs",
                                  mito_masks=mito_masks, lyso_masks=lyso_masks)
```

