
from skimage.measure import regionprops, label

def classify_mitochondria(binary_mask):
    """
    Classify mitochondria as elongated or punctate based on area and eccentricity.
    Updated stricter criteria:
    - Elongated: area >= 0.15 and eccentricity >= 0.85
    - Punctate: area >= 0.05 and eccentricity <= 0.70
    """
    labels = label(binary_mask > 0)
    props = regionprops(labels)

    elongated, punctate = [], []
    for prop in props:
        area = prop.area
        ecc = prop.eccentricity

        if area >= 0.15 and ecc >= 0.85:
            elongated.append(prop)
        elif area >= 0.05 and ecc <= 0.70:
            punctate.append(prop)

    return elongated, punctate
