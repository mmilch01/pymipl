#!/usr/bin/env python3
import argparse
import os

import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
from scipy.ndimage import center_of_mass
from nibabel.processing import resample_from_to


def load_and_canonical(path):
    """
    Load NIfTI image and reorient to RAS+ (orientation-proof step).
    This ensures all images/masks share a consistent orientation basis.
    """
    img = nib.load(path)
    img = nib.as_closest_canonical(img)
    return img


def resample_mask_to_image(mask_img, ref_img):
    """
    Resample mask to the reference image grid using nearest-neighbor
    interpolation to preserve binary labels and ensure proper
    alignment in both space and orientation.
    """
    if mask_img.shape != ref_img.shape or not np.allclose(mask_img.affine, ref_img.affine):
        mask_img = resample_from_to(mask_img, ref_img, order=0)  # nearest neighbor
    return mask_img


def compute_center_of_mass(mask_data):
    """
    Compute center of mass in voxel coordinates from a binary mask.
    Returns integer indices.
    """
    if np.count_nonzero(mask_data) == 0:
        raise ValueError("Mask 2 is empty; cannot compute center of mass.")
    com = center_of_mass(mask_data > 0)
    # center_of_mass returns floats; convert to nearest integer indices
    com_idx = tuple(int(round(c)) for c in com)
    return com_idx  # (i, j, k)


def get_three_plane_slices(image_data, mask1_data, mask2_data, com_idx):
    """
    Extract sagittal, coronal, and axial slices for image and masks
    through the center of mass.

    image_data, mask1_data, mask2_data: 3D arrays in the same space.
    com_idx: (i, j, k) voxel indices.
    """
    i, j, k = com_idx
    nx, ny, nz = image_data.shape

    # Clamp indices to valid range, just in case rounding pushed them out
    i = max(0, min(nx - 1, i))
    j = max(0, min(ny - 1, j))
    k = max(0, min(nz - 1, k))

    # Sagittal: slice along x (i fixed), plane is (y, z)
    sag_img = image_data[i, :, :]
    sag_m1 = mask1_data[i, :, :]
    sag_m2 = mask2_data[i, :, :]

    # Coronal: slice along y (j fixed), plane is (x, z)
    cor_img = image_data[:, j, :]
    cor_m1 = mask1_data[:, j, :]
    cor_m2 = mask2_data[:, j, :]

    # Axial: slice along z (k fixed), plane is (x, y)
    axi_img = image_data[:, :, k]
    axi_m1 = mask1_data[:, :, k]
    axi_m2 = mask2_data[:, :, k]

    return (sag_img, cor_img, axi_img), (sag_m1, cor_m1, axi_m1), (sag_m2, cor_m2, axi_m2)


def plot_three_plane(image_slices, mask1_slices, mask2_slices, output_png):
    """
    Create a PNG with three-plane views and overlay mask1 (green)
    and mask2 (red) as contours.
    """
    sag_img, cor_img, axi_img = image_slices
    sag_m1, cor_m1, axi_m1 = mask1_slices
    sag_m2, cor_m2, axi_m2 = mask2_slices

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    titles = ["Sagittal", "Coronal", "Axial"]

    # To keep orientation consistent between image and masks, we apply
    # the same rotation to each slice. Rot90 (k=1) is often convenient
    # for display; it does not break alignment as long as it is applied
    # identically to image and masks.
    data = [
        (sag_img, sag_m1, sag_m2),
        (cor_img, cor_m1, cor_m2),
        (axi_img, axi_m1, axi_m2),
    ]

    for ax, (img_slice, m1_slice, m2_slice), title in zip(axes, data, titles):
        img_disp = np.rot90(img_slice)
        m1_disp = np.rot90(m1_slice)
        m2_disp = np.rot90(m2_slice)

        ax.imshow(img_disp, cmap="gray", interpolation="nearest")

        # Contours: level 0.5 for binary masks (0/1)
        if np.any(m1_disp > 0):
            ax.contour(m1_disp, levels=[0.5], colors="g", linewidths=1)
        if np.any(m2_disp > 0):
            ax.contour(m2_disp, levels=[0.5], colors="r", linewidths=1)

        ax.set_title(title)
        ax.axis("off")

    plt.tight_layout()
    fig.savefig(output_png, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Create three-plane view of a NIfTI image with two binary mask contours, "
                    "centered on the center of mass of mask2. Uses orientation-proof alignment."
    )
    parser.add_argument(
        "image",
        type=str,
        help="Path to the NIfTI image (image 1)",
    )
    parser.add_argument(
        "mask1",
        type=str,
        help="Path to the first NIfTI binary mask (mask 1)",
    )
    parser.add_argument(
        "mask2",
        type=str,
        help="Path to the second NIfTI binary mask (mask 2, used for center of mass)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="three_plane_view.png",
        help="Output PNG filename (default: three_plane_view.png)",
    )

    args = parser.parse_args()

    # Load and canonicalize all images (orientation-proof step)
    img1 = load_and_canonical(args.image)
    mask1 = load_and_canonical(args.mask1)
    mask2 = load_and_canonical(args.mask2)

    # Resample masks to image grid/affine
    mask1 = resample_mask_to_image(mask1, img1)
    mask2 = resample_mask_to_image(mask2, img1)

    img_data = img1.get_fdata()
    mask1_data = mask1.get_fdata()
    mask2_data = mask2.get_fdata()

    if img_data.ndim != 3:
        raise ValueError("Image 1 must be 3D.")
    if mask1_data.shape != img_data.shape or mask2_data.shape != img_data.shape:
        raise ValueError("After resampling, masks must have the same shape as image.")

    # Ensure binary masks (in case they are not exactly 0/1)
    mask1_data = (mask1_data > 0).astype(np.uint8)
    mask2_data = (mask2_data > 0).astype(np.uint8)

    # Center of mass of mask2 in voxel indices
    com_idx = compute_center_of_mass(mask2_data)

    # Get three-plane slices
    image_slices, mask1_slices, mask2_slices = get_three_plane_slices(
        img_data, mask1_data, mask2_data, com_idx
    )

    # Generate figure
    output_png = args.output
    # Ensure directory exists if user provided a path
    out_dir = os.path.dirname(output_png)
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    plot_three_plane(image_slices, mask1_slices, mask2_slices, output_png)


if __name__ == "__main__":
    main()
