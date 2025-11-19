#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import numpy as np
import SimpleITK as sitk
from scipy.ndimage import distance_transform_edt


def load_mask(path):
    img = sitk.ReadImage(str(path))
    arr = sitk.GetArrayFromImage(img).astype(np.uint8)
    return arr


def dice_coefficient(a, b):
    a = a.astype(bool)
    b = b.astype(bool)
    inter = np.logical_and(a, b).sum()
    denom = a.sum() + b.sum()
    if denom == 0:
        return 1.0
    return 2 * inter / denom


def average_surface_distance(mask_a, mask_b):
    """
    Average symmetric surface distance:
    1. Extract surfaces of both binary masks.
    2. Compute EDT of each mask.
    3. Compute distances from surface A to mask B and from B to A.
    4. Return average.
    """

    # ensure boolean
    A = mask_a.astype(bool)
    B = mask_b.astype(bool)

    # surface = voxels on boundary: mask AND not eroded(mask)
    # (Simple and robust for binary masks)
    from scipy.ndimage import binary_erosion

    A_er = binary_erosion(A)
    B_er = binary_erosion(B)

    A_surf = np.logical_and(A, np.logical_not(A_er))
    B_surf = np.logical_and(B, np.logical_not(B_er))

    # distance transform of the opposite mask's interior
    dt_A = distance_transform_edt(~A)  # distance to nearest A
    dt_B = distance_transform_edt(~B)  # distance to nearest B

    # distances
    dist_A_to_B = dt_B[A_surf]
    dist_B_to_A = dt_A[B_surf]

    if dist_A_to_B.size == 0 and dist_B_to_A.size == 0:
        return 0.0

    vals = np.concatenate([dist_A_to_B, dist_B_to_A])
    return float(vals.mean())


def main():
    p = argparse.ArgumentParser(description="Compute Dice and average surface distance between two binary masks.")
    p.add_argument("--mask1", required=True, type=str, help="Path to first binary mask")
    p.add_argument("--mask2", required=True, type=str, help="Path to second binary mask")
    p.add_argument("--out", required=True, type=str, help="Output JSON file")
    
    args = p.parse_args()

    print (f"Comparing masks {args.mask1} and {args.mask2}")
    m1 = load_mask(args.mask1)
    m2 = load_mask(args.mask2)

    dice = dice_coefficient(m1, m2)
    asd = -1 #average_surface_distance(m1, m2)

    result = {
        "mask1": str(Path(args.mask1).resolve()),
        "mask2": str(Path(args.mask2).resolve()),
        "dice": float(dice),        
    }
    #"average_surface_distance": float(asd),

    with open(args.out, "w") as f:
        json.dump(result, f, indent=2)


if __name__ == "__main__":
    main()
