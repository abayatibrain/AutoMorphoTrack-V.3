
# AutoMorphoTrack

**AutoMorphoTrack** is a Python package for automated analysis of mitochondrial morphology, lysosomal count, organelle tracking, and colocalization from two-channel time-lapse microscopy images. It enables high-content, publication-ready quantification and visualization of organelle dynamics in live-cell or fixed imaging datasets.

---

## 📦 Features

- ✅ **Organelle Detection** (Lysosomes & Mitochondria)
- ✅ **Mitochondrial Morphology Classification** (Elongated vs. Punctate)
- ✅ **Lysosomal Count per Frame**
- ✅ **Organelle Tracking** (Displacement, Velocity)
- ✅ **Colocalization Analysis** (Mito-Lyso overlap using masked Pearson correlation)
- ✅ **Publication-Ready Image Output** (Annotated overlays and tracking maps)
- ✅ **Customizable Thresholds and Criteria**
- ✅ **Batch Compatibility for High-Throughput Analysis**

---

## 🧰 Installation

```bash
git clone https://github.com/yourusername/AutoMorphoTrack.git
cd AutoMorphoTrack
pip install -e .
```

This installs the package in **editable mode**, allowing you to edit `.py` files while keeping the package functional.

Dependencies are specified in the individual scripts and include:
- `numpy`
- `scipy`
- `opencv-python`
- `matplotlib`
- `scikit-image`
- `pandas`

You can install them via:

```bash
pip install numpy scipy opencv-python matplotlib scikit-image pandas
```

---

## 📂 Folder Structure

```
AutoMorphoTrack/
│
├── core.py                      # Pipeline controller
├── detection.py                # Organelle detection and segmentation
├── morphology.py               # Mitochondrial morphology classification
├── Lyso_count.py               # Lysosomal counting per frame
├── tracking.py                 # Organelle displacement & velocity
├── visualization.py            # All figure generation
├── AutoMorphoTrack_Full_Package.ipynb  # Primary usage notebook
├── AutoMorphoTrack_Example.ipynb       # Simplified example notebook
├── AutoMorphoTrack_Source_Code.ipynb   # Full code + explanations
├── README_AutoMorphoTrack.md   # This file
└── setup.py                    # Package installation script
```

---

## 🖼️ Input Format

- **Input**: 2-channel `.tif` time-lapse image stack
  - Channel 0: Lysosomes
  - Channel 1: Mitochondria
- **Dimensions**: (frames, height, width, channels)

---

## 🚀 How to Use

### 📓 Option 1: Notebook

Use the notebook `AutoMorphoTrack_Full_Package.ipynb` for interactive analysis.

```python
# Load and run this notebook cell by cell.
```

### 🧪 Option 2: Python Script

You can also build a script like below:

```python
from AutoMorphoTrack import core

core.run_pipeline(
    filepath="path/to/your/stack.tif",
    output_dir="path/to/results",
    lyso_channel=0,
    mito_channel=1,
    frame_to_visualize=0
)
```

---

## 🧮 Analysis Workflow

1. **Preprocessing**:
   - Adaptive thresholding for each channel
   - Watershed-based segmentation in dense regions

2. **Detection**:
   - Lysosomes: Detected as small, circular puncta
   - Mitochondria: Detected via area and shape features

3. **Morphology Classification**:
   - **Elongated**: area ≥ 0.2 µm², eccentricity ≥ 0.95
   - **Punctate**: area ≥ 0.015 µm², eccentricity ≤ 0.85

4. **Tracking**:
   - Organelle centroids linked across frames
   - Displacement and velocity calculated

5. **Colocalization**:
   - Masked Pearson correlation across overlapping signal regions
   - Output: scatterplot + colocalization overlay image

---

## 📊 Output Files

All output files are saved to the specified output directory:

### 1. **Images**
| File | Description |
|------|-------------|
| `threshold_mito_frame0.png` | Thresholded mitochondria |
| `threshold_lyso_frame0.png` | Thresholded lysosomes |
| `detection_frame0.png` | Raw image with detected organelles (mito = red, lyso = green) |
| `morphology_frame0.png` | Labels each mitochondrion with `E` (elongated) or `P` (punctate) |
| `lysosome_labels_frame0.png` | Each lysosome labeled with a number |
| `colocalization_frame0.png` | Mitochondria and lysosomes shown with cyan overlap |
| `tracks_mito_frame0.png` | All mitochondrial tracks overlaid on frame 0 |
| `tracks_lyso_frame0.png` | All lysosomal tracks overlaid on frame 0 |

### 2. **CSV Tables**
| File | Description |
|------|-------------|
| `morphology_summary.csv` | Elongated vs punctate per frame |
| `lysosome_count.csv` | Lysosomal count per frame |
| `displacement_velocity_mito.csv` | Mitochondrial tracking statistics |
| `displacement_velocity_lyso.csv` | Lysosomal tracking statistics |
| `pearson_colocalization.csv` | Pearson r value per frame (masked) |

---

## 🎛️ Customization

You can modify the following in `core.py` or the notebooks:

- **Thresholding method** (adaptive, global, or Otsu)
- **Segmentation settings** (`min_size`, `circularity`, etc.)
- **Morphology criteria** (eccentricity/area thresholds)
- **Frame of interest** for visualization

---

## 📑 Manuscript Integration

This package is designed for inclusion in bioimage analysis pipelines in publications. Figures generated here are:

- Resolution-independent (via Matplotlib)
- Consistent across frames and conditions
- Automatically saved for direct inclusion in manuscripts

---

## 🧪 Testing

Run the test script to verify basic functionality:

```bash
python test_basic.py
```

---

## 🧠 Citation

If you use AutoMorphoTrack in your research, please cite:

**Bayati, A.**, *et al.* (2025). *AutoMorphoTrack: A high-throughput pipeline for tracking and classification of lysosomes and mitochondria from time-lapse microscopy*. (Preprint/DOI TBD)

---

## 📬 Contact

For questions, issues, or contributions, contact:  
**Armin Bayati**  
Research Fellow, Harvard Medical School / MGH  
[GitHub Profile](https://github.com/arminbayati)


## 📊 Additional Quantification Modules
- `shape_descriptors.csv`: Mitochondrial area, eccentricity, circularity, solidity, aspect ratio
- `aspect_ratio_boxplot.png`: Boxplot of aspect ratio by frame
- `circularity_boxplot.png`: Boxplot of circularity by frame
- `directionality_summary.csv`: Directionality index per mitochondrion
- `directionality_hist.png`: Histogram of directionality index
- `lyso_mito_min_distances.csv`: Distance from each lysosome to closest mitochondrion
- `lyso_mito_distance_hist.png`: Histogram of minimum lyso-mito distances


## 🔄 Default Detection Settings
- As of the latest version, AutoMorphoTrack uses **relaxed detection criteria**:
  - `min_size = 5` pixels for both mitochondria and lysosomes
  - This allows for better sensitivity, especially in densely populated fields


## 📐 Organelle Shape Descriptor Analysis
- Calculates circularity and aspect ratio across all frames
- Generates CSVs and visualizations with trendlines
- Helps detect fragmentation, elongation, or rounding trends


### 🔍 Improved Organelle Detection Overlay
- Shows grayscale images with red outlines of detected structures
- Available for both mitochondria and lysosomes


### July 2025 Update
- Updated mitochondrial morphology classification:
  - Elongated: area ≥ 0.1 µm², eccentricity ≥ 0.7
  - Punctate: area ≥ 0.01 µm², eccentricity ≤ 0.4
- Fixed a bug in visualization that attempted to render extremely large image sizes.
- Added color-coded dot visualization for organelle tracking overlayed on Frame 0.
  - Output supported for both mitochondria and lysosomes.
- Added automatic export of circularity and aspect ratio metrics as CSV.
- Added default graphing of these shape descriptors across time as PNG plots.
- Added trendline overlay to circularity and aspect ratio plots.
- Shape analysis now automatically runs for both mitochondria and lysosomes.
- Added masked Pearson correlation across frames as default colocalization output.
- Added visualization of colocalized (cyan) pixels on 5 representative frames.
- Colocalization overlays are now compiled into a montage for easier visualization.
- Organelle tracking module now includes 4 visualization options:
  - Arrowed trajectories
  - Rainbow gradient trails
  - Heatmap of organelle occupancy
  - Start/end point summary
- All visualizations are automatically saved as PNG files.
- Morphology classification now exports per-frame elongated and punctate counts to CSV and PNG.

## AutoMorphoTrack - Full Pipeline

This package includes a complete analysis pipeline for mitochondrial and lysosomal tracking, morphology classification, and motion analysis. The main notebook `AutoMorphoTrack_Full_Package_Updated.ipynb` runs the following steps:

1. Load a multi-channel .tif image stack
2. Assign channel 0 as mitochondria and channel 1 as lysosomes
3. Apply threshold-based detection for both organelles
4. Classify mitochondrial morphology (elongated vs punctate)
5. Compute displacement and velocity across frames
6. Output:
   - CSVs: `mitochondria_morphology.csv`, `mitochondria_motion.csv`, `lysosome_motion.csv`
   - PNGs: (optional if added for plots)


## New Feature: Full Pipeline
You can now run the entire analysis using:
```python
from analysis import run_full_analysis
run_full_analysis(image_path, output_dir, mito_channel=0, lyso_channel=1)
```


### Updated Pipeline
The analysis now uses `detection.py` for organelle detection (Otsu threshold + labeling).
This improves segmentation and enables more accurate morphology and motion analysis.


### Detection Output Update
- Now also saves individual grayscale channel images:
  - `Channel_0_Lysosomes_Frame0.png`
  - `Channel_1_Mitochondria_Frame0.png`


### Morphology Criteria Update
- Elongated: area ≥ 0.02, eccentricity ≥ 0.75
- Punctate: area ≥ 0.025, eccentricity ≤ 0.65
- Font size for morphology labels increased to 10


### Lysosome Counting Module
- Counts lysosomes per frame using Otsu thresholding and connected components
- Outputs: count CSV, labeled images, and a count plot
