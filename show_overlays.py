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
    com_idx = tuple(int(round(c)) for c in com)
    return com_idx  # (i, j, k)


def get_three_plane_slices(image_data, mask1_data, mask2_data, com_idx, mask3_data=None):
    """
    Extract sagittal, coronal, and axial slices for image and masks
    through the center of mass.

    image_data, mask1_data, mask2_data, mask3_data: 3D arrays in the same space.
    com_idx: (i, j, k) voxel indices.
    """
    i, j, k = com_idx
    nx, ny, nz = image_data.shape

    # Clamp indices
    i = max(0, min(nx - 1, i))
    j = max(0, min(ny - 1, j))
    k = max(0, min(nz - 1, k))

    # Sagittal: slice along x (i fixed), plane is (y, z)
    sag_img = image_data[i, :, :]
    sag_m1 = mask1_data[i, :, :]
    sag_m2 = mask2_data[i, :, :]
    sag_m3 = mask3_data[i, :, :] if mask3_data is not None else None

    # Coronal: slice along y (j fixed), plane is (x, z)
    cor_img = image_data[:, j, :]
    cor_m1 = mask1_data[:, j, :]
    cor_m2 = mask2_data[:, j, :]
    cor_m3 = mask3_data[:, j, :] if mask3_data is not None else None

    # Axial: slice along z (k fixed), plane is (x, y)
    axi_img = image_data[:, :, k]
    axi_m1 = mask1_data[:, :, k]
    axi_m2 = mask2_data[:, :, k]
    axi_m3 = mask3_data[:, :, k] if mask3_data is not None else None

    return (
        (sag_img, cor_img, axi_img),
        (sag_m1, cor_m1, axi_m1),
        (sag_m2, cor_m2, axi_m2),
        (sag_m3, cor_m3, axi_m3),
    )

def plot_three_plane(image_slices, mask1_slices, mask2_slices, output_png,
                     mask3_slices=None, spacing=None, palette="green,red,lightblue"):
    """
    Create a PNG with three-plane views and overlay mask1 (green),
    mask2 (red), and optional mask3 (blue) as contours.
    Aspect ratio is preserved using voxel spacing from input image.
    """

    if spacing is None:
        raise ValueError("Voxel spacing array must be provided.")
    dx, dy, dz = spacing

    # unpack slices
    sag_img, cor_img, axi_img = image_slices
    sag_m1, cor_m1, axi_m1 = mask1_slices
    sag_m2, cor_m2, axi_m2 = mask2_slices

    if mask3_slices is not None:
        sag_m3, cor_m3, axi_m3 = mask3_slices
    else:
        sag_m3 = cor_m3 = axi_m3 = None

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    titles = ["Sagittal", "Coronal", "Axial"]

    data = [

        # SAGITTAL slice: data is (Y,Z) = dy,dz spacing
        (
            np.rot90(sag_img),
            np.rot90(sag_m1),
            np.rot90(sag_m2),
            np.rot90(sag_m3) if sag_m3 is not None else None,
            [0, sag_img.shape[0] * dy, 0, sag_img.shape[1] * dz],
        ),

        # CORONAL slice: data is (X,Z) = dx,dz spacing
        (
            np.rot90(cor_img),
            np.rot90(cor_m1),
            np.rot90(cor_m2),
            np.rot90(cor_m3) if cor_m3 is not None else None,
            [0, cor_img.shape[0] * dx, 0, cor_img.shape[1] * dz],
        ),

        # AXIAL slice: data is (X,Y) = dx,dy spacing
        (
            np.rot90(axi_img),
            np.rot90(axi_m1),
            np.rot90(axi_m2),
            np.rot90(axi_m3) if axi_m3 is not None else None,
            [0, axi_img.shape[0] * dx, 0, axi_img.shape[1] * dy],
        ),
    ]

    colors=palette.split(',')
    
    for ax, (img_disp, m1_disp, m2_disp, m3_disp, extent), title in zip(axes, data, titles):

        ax.imshow(img_disp, cmap="gray", interpolation="nearest",
                  extent=extent, aspect='equal',origin='upper')

        if np.any(m1_disp > 0):
            ax.contour(m1_disp, levels=[0.5], colors=colors[0],
                       linewidths=1, extent=extent,origin='upper')

        if np.any(m2_disp > 0):
            ax.contour(m2_disp, levels=[0.5], colors=colors[1],
                       linewidths=1, extent=extent,origin='upper')

        if m3_disp is not None and np.any(m3_disp > 0):
            ax.contour(m3_disp, levels=[0.5], colors=colors[2],
                       linewidths=1, extent=extent,origin='upper')

        #ax.set_title(title)
        ax.axis("off")

    plt.tight_layout()
    fig.savefig(output_png, dpi=150, bbox_inches="tight")
    plt.close(fig)



def main():
    parser = argparse.ArgumentParser(
        description=(
            "Create three-plane view of a NIfTI image with two or three binary mask contours, "
            "centered on the center of mass of mask2. Uses orientation-proof alignment."
        )
    )
    parser.add_argument(
        "image",
        type=str,
        help="Path to the NIfTI image (image 1)",
    )
    parser.add_argument(
        "mask1",
        type=str,
        help="Path to the first NIfTI binary mask (mask 1, green)",
    )
    parser.add_argument(
        "mask2",
        type=str,
        help="Path to the second NIfTI binary mask (mask 2, red; used for center of mass)",
    )
    parser.add_argument(
        "--mask3",
        type=str,
        default=None,
        help="Optional third NIfTI binary mask (mask 3, light blue)",
    )
    parser.add_argument(
        "--palette",
        type=str,
        default="green,red,lighblue",
        help="mask drawing palette, comma separated [green,red,lightblue]"
    )
    
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="three_plane_view.png",
        help="Output PNG filename (default: three_plane_view.png)",
    )
    parser.add_argument(
        "--max-dim",
        type=int,
        default=256,
        help="Resize image/masks so largest dimension <= max-dim (default 256).",
    )

    args = parser.parse_args()

    
    # Load and canonicalize main image and required masks
    img1 = load_and_canonical(args.image)
    mask1 = load_and_canonical(args.mask1)
    mask2 = load_and_canonical(args.mask2)

    # Optional mask3
    mask3 = None
    print (f"{__file__}: loading {args.image}, {args.mask1}, {args.mask2}")
    
    if args.mask3 is not None:
        print (f"{__file__}: loading {args.mask3}")
        mask3 = load_and_canonical(args.mask3)
        print ('Computing three-mask overlay')
    else:
        print ('Computing two-mask overlay')

    # Resample masks to image grid/affine
    mask1 = resample_mask_to_image(mask1, img1)
    mask2 = resample_mask_to_image(mask2, img1)
    if mask3 is not None:
        mask3 = resample_mask_to_image(mask3, img1)

    img_data = img1.get_fdata()
    mask1_data = mask1.get_fdata()
    mask2_data = mask2.get_fdata()
    mask3_data = mask3.get_fdata() if mask3 is not None else None

    if img_data.ndim != 3:
        raise ValueError("Image 1 must be 3D.")
    if mask1_data.shape != img_data.shape or mask2_data.shape != img_data.shape:
        raise ValueError("After resampling, mask1 and mask2 must have the same shape as image.")
    if mask3_data is not None and mask3_data.shape != img_data.shape:
        raise ValueError("After resampling, mask3 must have the same shape as image.")

    # Ensure binary masks
    mask1_data = (mask1_data > 0).astype(np.uint8)
    mask2_data = (mask2_data > 0).astype(np.uint8)
    if mask3_data is not None:
        mask3_data = (mask3_data > 0).astype(np.uint8)
    
    # Center of mass of mask2
    com_idx = compute_center_of_mass(mask2_data)

    # Get three-plane slices (image + masks)
    image_slices, mask1_slices, mask2_slices, mask3_slices = get_three_plane_slices(
        img_data, mask1_data, mask2_data, com_idx, mask3_data=mask3_data
    )

    output_png = args.output
    out_dir = os.path.dirname(output_png)
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    plot_three_plane(
        image_slices,
        mask1_slices,
        mask2_slices,
        output_png,
        mask3_slices=mask3_slices if mask3_data is not None else None,
        spacing=img1.header.get_zooms()[:3],
        palette=args.palette
    )


if __name__ == "__main__":
    main()
