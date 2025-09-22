
# AutoMorphoTrack

**AutoMorphoTrack** is a fully automated pipeline for analyzing mitochondrial and lysosomal dynamics, morphology, and interactions from multi-channel time-lapse fluorescence microscopy stacks.

---

## 🔧 Features and Capabilities

- ✅ Automatic detection and thresholding for mitochondria (Channel 0) and lysosomes (Channel 1)
- 🔍 Morphology classification: punctate vs. elongated mitochondria (with visual overlays and CSV export)
- 📈 Motion analysis: displacement and velocity calculated per organelle, per frame
- 🔄 Vector visualization of organelle motion (displacement arrows on Frame 0)
- 🧠 Shape feature extraction: area, eccentricity, circularity, solidity, orientation, aspect ratio
- 📊 CSV output + automated charting of:
  - Morphology counts
  - Displacement and velocity
  - Shape features
- 🔗 Colocalization analysis using masked Pearson correlation

---

## 📁 Input Requirements

- Multi-channel `.tif` stack
- Channel 0: Mitochondria  
- Channel 1: Lysosomes  
- Channel 2: (ignored if present)

---

## 🚀 How to Use (Python Users)

```python
from automorphotrack import run_full_pipeline

run_full_pipeline("your_image_stack.tif")
```

Outputs are saved automatically to `./outputs/` including figures, overlays, CSVs, and charts.

---

## 🤖 Using with AI Chatbots (Non-Coding Users)

Ask an AI assistant like ChatGPT to:

> “Run the AutoMorphoTrack pipeline on `example.tif`, mitochondria in channel 0, lysosomes in channel 1.”

It will walk you through detection, classification, motion analysis, and generate publication-ready outputs.

---

## 📦 Outputs

| Output Type         | Description                                   |
|---------------------|-----------------------------------------------|
| Figures             | Overlay images, vector motion, shape plots    |
| CSV Files           | Morphology counts, shape features, velocities |
| Charts              | Line plots for shape, motion, colocalization  |

---

## 🧪 Citation

If you use AutoMorphoTrack in your work, please cite:
> Bayati, A. et al. AutoMorphoTrack: automated analysis of organelle dynamics and morphology. *In preparation*.

---

## 🖼️ Colocalization Visualization (NEW)

AutoMorphoTrack now includes **visualization of mitochondria–lysosome colocalization** in addition to masked Pearson correlation values.

### What you get
- **Frame 0 composite + colocalization heatmap** (`colocalization_frame0.png`):  
  - Left: RGB composite (R = mitochondria, G = lysosomes)  
  - Right: Colocalization heatmap computed as normalized intensity product
- **Time series plot** (`colocalization_timeseries.png`): masked Pearson *r* per frame
- **CSV output** (`colocalization_timeseries.csv`): frame-wise masked Pearson *r*

### Minimal example

```python
import tifffile as tiff
import numpy as np
from colocalization import run_colocalization_pipeline

# Load your two-channel stack (T x H x W x C) or (C x T x H x W); adapt as needed
stack = tiff.imread("path/to/stack.tif")

# Example: if stack is (T, H, W, 2)
mito_stack = stack[..., 0]
lyso_stack = stack[..., 1]

out = run_colocalization_pipeline(mito_stack, lyso_stack, output_dir="outputs")
print(out)
```
Outputs are saved under `outputs/`:
- `colocalization_frame0.png`
- `colocalization_timeseries.png`
- `colocalization_timeseries.csv`

> **Note:** If you already compute masks for mitochondria/lysosomes, you can pass them as `mito_masks` and `lyso_masks` to restrict correlation to organelle-positive pixels.



### Lysosomal Count Analysis

AutoMorphoTrack now automatically performs lysosomal counting:

- **`lysosome_counts.csv`** → frame-by-frame lysosomal counts  
- **`lysosome_count_overlay_Frame0.png`** → annotated Frame 0 image with count labels  
- **`lysosome_count_over_time.png`** → line graph of lysosome count across time  

This analysis runs **immediately after mitochondrial morphology classification** in the standard pipeline.


### Mitochondrial Morphology Classification (Updated)

Elongated mitochondria are now classified with **relaxed thresholds** to increase detection:
- **Elongated:** area ≥ 0.15 µm² and eccentricity ≥ 0.7  
- **Punctate:** area ≥ 0.015 µm² and eccentricity ≤ 0.85  

This change increases the number of mitochondria detected as elongated while keeping punctate classification stable.
