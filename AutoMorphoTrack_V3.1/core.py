from skimage.segmentation import watershed
from scipy import ndimage as ndi
from skimage.morphology import dilation, disk

from .detection import detect_organelles
from .tracking import calculate_displacement
from .morphology import classify_mitochondria
from .visualization import plot_organelles

def run_full_pipeline(image_stack):
    for frame in image_stack:
        mask = detect_organelles(frame)
        plot_organelles(frame, mask)



import tifffile
import numpy as np

def load_stack(path):
    """
    Load a TIFF stack and standardize to shape (frames, channels, height, width).
    Supports TIFFs with channels as last dimension or second dimension.
    """
    stack = tifffile.imread(path)

    if stack.ndim == 4:
        if stack.shape[-1] <= 4:  # format (frames, H, W, C)
            stack = np.transpose(stack, (0, 3, 1, 2))
        elif stack.shape[1] <= 4:  # format (frames, C, H, W)
            pass  # already correct
        else:
            raise ValueError(f"Unexpected stack shape: {stack.shape}")
    else:
        raise ValueError(f"Unsupported stack dimensions: {stack.shape}")

    return stack
