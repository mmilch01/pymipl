#!/usr/bin/env python3
import argparse
import os

import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
from scipy.ndimage import center_of_mass
from nibabel.processing import resample_from_to

PALETTE = "green,red,lightblue,yellow,magenta,cyan,orange,purple,brown,pink,olive,teal,navy,gold,gray,black"

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

def get_three_plane_slices(image_data, mask_data, com_idx=None):
    """
    Extract sagittal, coronal, and axial slices for image and a list of masks.

    image_data: 3D array
    masks_data: list[3D array] in the same space
    com_idx: (i, j, k) voxel indices; if None, uses central slices.
    """
    nx, ny, nz = image_data.shape

    if com_idx is None:
        i, j, k = nx // 2, ny // 2, nz // 2
    else:
        i, j, k = com_idx

    # Clamp indices
    i = max(0, min(nx - 1, i))
    j = max(0, min(ny - 1, j))
    k = max(0, min(nz - 1, k))

    # Image slices
    sag_img = image_data[i, :, :]
    cor_img = image_data[:, j, :]
    axi_img = image_data[:, :, k]
    image_slices = (sag_img, cor_img, axi_img)

    # Mask slices list
    masks_slices = []
    for m in mask_data:
        masks_slices.append((m[i, :, :], m[:, j, :], m[:, :, k]))

    return image_slices, masks_slices
    
def plot_three_plane(image_slices, masks_slices, output_png, spacing=None, palette=PALETTE):
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
    
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), gridspec_kw={"wspace": 0, "hspace": 0})

    titles = ["Sagittal", "Coronal", "Axial"]
    
    # pre-rotate all mask slices per plane (same rotation as image)
    sag_masks = [np.rot90(m[0]) for m in masks_slices]
    cor_masks = [np.rot90(m[1]) for m in masks_slices]
    axi_masks = [np.rot90(m[2]) for m in masks_slices]
    
    data = [
        # SAGITTAL (Y,Z)
        (np.rot90(sag_img), sag_masks, [0, sag_img.shape[0] * dy, 0, sag_img.shape[1] * dz]),
        # CORONAL (X,Z)
        (np.rot90(cor_img), cor_masks, [0, cor_img.shape[0] * dx, 0, cor_img.shape[1] * dz]),
        # AXIAL (X,Y)
        (np.rot90(axi_img), axi_masks, [0, axi_img.shape[0] * dx, 0, axi_img.shape[1] * dy]),
    ]
    
    colors = [c.strip() for c in palette.split(",") if c.strip()]

    #compute robust min and max
    all_vals = np.concatenate([sag_img.ravel(), cor_img.ravel(),   axi_img.ravel() ])
    all_vals = all_vals[np.isfinite(all_vals) & (all_vals != 0)]
    if all_vals.size == 0:
        #raise ValueError("Base image has no finite non-zero voxels.")
        vmin, vmax = 0,0
    vmin, vmax = np.percentile(all_vals, [2, 98])
    
    for ax, (img_disp, masks_disp, extent), title in zip(axes, data, titles):
        ax.imshow(img_disp, cmap="gray", interpolation="nearest",
                  extent=extent, aspect='equal', origin='upper',vmin=vmin,vmax=vmax)
    
        for idx, m in enumerate(masks_disp):
            if np.any(m > 0) and colors:
                ax.contour(m, levels=[0.5], colors=colors[idx % len(colors)],
                           linewidths=1, extent=extent, origin='upper')

        ax.axis("off")

    fig.subplots_adjust(left=0, right=1, top=1, bottom=0, wspace=0, hspace=0)
    fig.savefig(output_png, dpi=150, bbox_inches="tight", pad_inches=0)
    plt.close(fig)



def main():
    parser = argparse.ArgumentParser(
        description=(
            "Create three-plane view of a NIfTI image with binary mask contours, "
            "centered on the center of mass of mask 1 (if specified). Uses orientation-proof alignment."
        )
    )
    parser.add_argument(
        "image",
        type=str,
        help="Path to the NIfTI image (image 1)",
    )
    parser.add_argument(
        "--mask",
        action="append",
        default=[],
        help="Repeatable NIfTI binary mask path. Use multiple times: --mask a.nii.gz --mask b.nii.gz ..."
    )
    parser.add_argument(
        "--palette",
        type=str,
        default="green,red,lighblue,...",
        help="mask drawing palette, comma separated (max 16 colors)"
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

    img1 = load_and_canonical(args.image)
    
    masks = []
    print(f"{__file__}: loading {args.image}")
    for p in args.mask:
        print(f"{__file__}: loading {p}")
        masks.append(load_and_canonical(p))
    
    print(f"Computing overlay with {len(masks)} mask(s)")
    

    # Resample masks to image grid/affine
    masks = [resample_mask_to_image(m, img1) for m in masks]

    img_data = img1.get_fdata()
    masks_data = [m.get_fdata() for m in masks]

    if img_data.ndim != 3:
        raise ValueError("Image 1 must be 3D.")

    for idx, m in enumerate(masks_data):
        if m.shape != img_data.shape:
            raise ValueError(f"After resampling, mask #{idx+1} must have the same shape as image.")
        
    # Ensure binary masks
    masks_data = [(m > 0).astype(np.uint8) for m in masks_data]
    
    # Center of mass of first mask (if any); otherwise use central slices
    com_idx = compute_center_of_mass(masks_data[0]) if masks_data else None

    image_slices, masks_slices = get_three_plane_slices(img_data, masks_data, com_idx=com_idx)

    output_png = args.output
    out_dir = os.path.dirname(output_png)
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    plot_three_plane(
        image_slices,
        masks_slices,
        output_png,
        spacing=img1.header.get_zooms()[:3],
        palette=args.palette
    )

if __name__ == "__main__":
    main()
