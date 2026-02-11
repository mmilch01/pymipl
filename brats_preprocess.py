#!/usr/bin/env python3
import os
import shutil
import subprocess
import argparse
import ants
import dicom2nifti

def brats_preprocess_pipeline(patient_id: str, modality_folders: dict, output_dir: str) -> None:
    """
    Convert DICOM->NIfTI (RAI), skull-strip reference (prefer T1) via HD-BET, rigid-register others to ref,
    apply ref brain mask, write outputs to output_dir/final/.
    """
    temp_dir = os.path.join(output_dir, "temp")
    final_dir = os.path.join(output_dir, "final")
    if os.path.isdir(temp_dir):
        shutil.rmtree(temp_dir)
    if os.path.isdir(final_dir):
        shutil.rmtree(final_dir)
    os.makedirs(temp_dir, exist_ok=True)
    os.makedirs(final_dir, exist_ok=True)

    # keep only provided modalities
    modality_folders = {k: v for k, v in modality_folders.items() if v}
    if not modality_folders:
        raise ValueError("No modality DICOM directories provided.")

    # 1) DICOM -> NIfTI + reorient to RAI
    nifti_files = {}
    for mod, dicom_path in modality_folders.items():
        out_file = os.path.join(temp_dir, f"{mod}.nii.gz")
        dicom2nifti.dicom_series_to_nifti(dicom_path, out_file, reorient_nifti=True)
        img = ants.image_read(out_file)
        img = ants.reorient_image2(img, "RAI")
        ants.image_write(img, out_file)
        nifti_files[mod] = out_file

    # 2) Choose reference (prefer T1)
    ref_mod = "t1" if "t1" in nifti_files else next(iter(nifti_files))
    ref_brain_path = os.path.join(final_dir, f"{patient_id}_{ref_mod}_brain.nii.gz")

    # 3) Skull strip reference using HD-BET
    subprocess.run(
        ["python", "-m", "HD_BET.entry_point", "--save_bet_mask", "-i", nifti_files[ref_mod], "-o", ref_brain_path],
        capture_output=True,
        text=True,
        check=True
    )
    ref_mask_path = ref_brain_path.replace(".nii.gz", "_bet.nii.gz")
    ref_mask = ants.image_read(ref_mask_path)
    ref_brain = ants.image_read(ref_brain_path)

    # 4) Register + mask other modalities to reference
    for mod, moving_path in nifti_files.items():
        if mod == ref_mod:
            continue
        moving_img = ants.image_read(moving_path)
        reg = ants.registration(fixed=ref_brain, moving=moving_img, type_of_transform="Rigid")
        stripped = ants.mask_image(reg["warpedmovout"], ref_mask)
        ants.image_write(stripped, os.path.join(final_dir, f"{patient_id}_{mod}_brain.nii.gz"))

    shutil.rmtree(temp_dir)


def _parse_args():
    p = argparse.ArgumentParser(description="BraTS-style preprocessing from explicit DICOM series directories.")
    p.add_argument("--patient_id", required=True)
    p.add_argument("--outdir", required=True, help="Output directory (creates temp/ and final/ inside).")
    p.add_argument("--t1", help="Path to T1 DICOM series directory.")
    p.add_argument("--t1ce", help="Path to T1 post-contrast (T1ce) DICOM series directory.")
    p.add_argument("--t2", help="Path to T2 DICOM series directory.")
    p.add_argument("--flair", help="Path to FLAIR DICOM series directory.")
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    mods = {"t1": args.t1, "t1ce": args.t1ce, "t2": args.t2, "flair": args.flair}
    brats_preprocess_pipeline(args.patient_id, mods, args.outdir)