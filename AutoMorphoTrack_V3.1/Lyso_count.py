import pandas as pd
from skimage.measure import label, regionprops

def count_lysosomes(mask_stack, output_dir):
    """
    Count lysosomes per frame and save results as CSV.
    """
    counts = []
    for frame_idx, mask in enumerate(mask_stack):
        labeled = label(mask)
        props = regionprops(labeled)
        counts.append({"Frame": frame_idx, "Lysosome_Count": len(props)})

    df = pd.DataFrame(counts)
    csv_path = f"{output_dir}/lysosome_counts.csv"
    df.to_csv(csv_path, index=False)
    return df


from skimage.measure import regionprops, label
import cv2

def detect_lysosomes(binary_mask, image_gray):
    """
    Detect lysosomes with stricter defaults and overlay labels.
    Criteria:
    - Area >= 8 pixels
    - Solidity >= 0.80
    Visualization:
    - Small font size (0.3), green numbering
    """
    labels = label(binary_mask > 0)
    props = regionprops(labels)

    filtered_props = [p for p in props if p.area >= 8 and p.solidity >= 0.80]

    # Create overlay
    overlay = cv2.cvtColor(image_gray, cv2.COLOR_GRAY2RGB)
    for i, prop in enumerate(filtered_props, start=1):
        y, x = prop.centroid
        cv2.putText(overlay, str(i), (int(x), int(y)), cv2.FONT_HERSHEY_SIMPLEX,
                    0.3, (0, 255, 0), 1, cv2.LINE_AA)

    return overlay, filtered_props
