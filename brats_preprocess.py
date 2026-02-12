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
    raw_dir = os.path.join(output_dir, "raw")
    bet_dir = os.path.join(output_dir, "bet")
    if os.path.isdir(raw_dir):
        shutil.rmtree(raw_dir)
    if os.path.isdir(bet_dir):
        shutil.rmtree(bet_dir)
    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(bet_dir, exist_ok=True)

    # keep only provided modalities
    modality_folders = {k: v for k, v in modality_folders.items() if v}
    if not modality_folders:
        raise ValueError("No modality DICOM directories provided.")

    # 1) DICOM -> NIfTI + reorient to RAI
    nifti_files = {}
    mod_spacing = {}
    for mod, dicom_path in modality_folders.items():
        print(f'converting DICOM to NIFTI for {mod}')
        out_file = os.path.join(raw_dir, f"{mod}.nii.gz")
        dicom2nifti.dicom_series_to_nifti(dicom_path, out_file, reorient_nifti=True)
        img = ants.image_read(out_file)
        img = ants.reorient_image2(img, "RAI")
        ants.image_write(img, out_file)
        nifti_files[mod] = out_file
        mod_spacing[mod] = img.spacing

    # 2) Choose reference with maximum minimum resolution 
    candidates = [m for m in ("t1", "t1ce") if m in nifti_files]
    if candidates:
        ref_mod = max(candidates, key=lambda m: min(mod_spacing[m]))
    else:
        ref_mod = next(iter(nifti_files))    
    ref_brain_path = os.path.join(bet_dir, f"{patient_id}_{ref_mod}_brain.nii.gz")

    # 3) Skull strip reference using HD-BET
    print('running HD-BET')
    subprocess.run(
        ["python", "-m", "HD_BET.entry_point", "--save_bet_mask", "-i", nifti_files[ref_mod], "-o", ref_brain_path],
        capture_output=True,
        text=True,
        check=True
    )
    ref_mask_path = ref_brain_path.replace(".nii.gz", "_bet.nii.gz")
    os.symlink(os.path.basename(ref_mask_path), os.path.join(os.path.dirname(os.path.abspath(ref_mask_path)),'brain_mask.nii.gz'))
               
    ref_mask = ants.image_read(ref_mask_path)
    ref_brain = ants.image_read(ref_brain_path)

    # 4) Register + mask other modalities to reference
    print('Running spatial registration and applying brain mask to other modalities')
    for mod, moving_path in nifti_files.items():
        if mod == ref_mod:
            continue
        moving_img = ants.image_read(moving_path)
        reg = ants.registration(fixed=ref_brain, moving=moving_img, type_of_transform="Rigid")
        stripped = ants.mask_image(reg["warpedmovout"], ref_mask)
        ants.image_write(stripped, os.path.join(bet_dir, f"{patient_id}_{mod}_brain.nii.gz"))
    
    #shutil.rmtree(raw_dir)


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