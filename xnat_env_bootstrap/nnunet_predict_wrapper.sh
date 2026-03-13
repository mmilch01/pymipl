#!/bin/bash
set -eu

usage() {
    cat <<'EOF'
Usage:
  nnunetv2_predict_wrapper.sh <input_nifti> <output_nifti> <dataset_id> [options]

Positional arguments:
  input_nifti         Input NIfTI file
  output_nifti        Output NIfTI file
  dataset_id          Dataset ID (folder name under $nnUNet_results)

Options:
  --nproc <n>         Number of processes for both -npp and -nps (default: 3)
  --use-gpu <0|1>     1 => device=cuda, 0 => device=cpu (default: 0)

Assumptions enforced:
  - $nnUNet_results/<dataset_id> contains exactly one configuration subfolder
  - configuration subfolder name is <trainer>__<plans>__<configuration>
  - all fold_* subfolders are used
  - each fold_* contains exactly one checkpoint_* file
  - all folds must have the same checkpoint filename
EOF
}

exit_with_error() {
    echo "ERROR: $*" >&2
    exit 1
}

[[ $# -ge 3 ]] || { usage; exit 1; }

input_file="$1"
output_file="$2"
dataset_id="$3"
shift 3

nproc=3
use_gpu=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --nproc) [[ $# -ge 2 ]] || exit_with_error "Missing value for $1"; nproc="$2"; shift 2 ;;
        --use-gpu) [[ $# -ge 2 ]] || exit_with_error "Missing value for $1"; use_gpu="$2"; shift 2 ;;
        *) exit_with_error "Unknown argument: $1" ;;
    esac
done

[[ -f "$input_file" ]] || exit_with_error "Input file does not exist: $input_file"
[[ "$use_gpu" == "0" || "$use_gpu" == "1" ]] || exit_with_error "--use-gpu must be 0 or 1"
[[ "$nproc" =~ ^[0-9]+$ ]] || exit_with_error "--nproc must be a non-negative integer"
[[ -n "${nnUNet_results}" ]] || exit_with_error "Environment variable nnUNet_results is not set"

dataset_dir="${nnUNet_results%/}/$dataset_id"
[[ -d "$dataset_dir" ]] || exit_with_error "Dataset directory does not exist: $dataset_dir"

mapfile -t config_dirs < <(find "$dataset_dir" -mindepth 1 -maxdepth 1 -type d | sort)
[[ "${#config_dirs[@]}" -eq 1 ]] || exit_with_error "Expected exactly one configuration folder in $dataset_dir, found ${#config_dirs[@]}"

config_dir="${config_dirs[0]}"
config_base="$(basename "$config_dir")"

trainer_id="${config_base%%__*}"; rest="${config_base#*__}"
plans_id="${rest%%__*}"
configuration_id="${rest#*__}"


[[ -n "${trainer_id}" ]] || exit_with_error "Could not parse trainer ID from configuration folder name: $config_base"
[[ -n "${plans_id}" ]] || exit_with_error "Could not parse plans ID from configuration folder name: $config_base"
[[ -n "${configuration_id}" ]] || exit_with_error "Could not parse configuration ID from configuration folder name: $config_base"


for required_file in dataset_fingerprint.json dataset.json plans.json; do
    [[ -f "$config_dir/$required_file" ]] || exit_with_error "Missing required file: $config_dir/$required_file"
done

mapfile -t fold_dirs < <(find "$config_dir" -mindepth 1 -maxdepth 1 -type d -name 'fold_*' | sort)
[[ "${#fold_dirs[@]}" -gt 0 ]] || exit_with_error "No fold_* folders found in $config_dir"

fold_ids=()
checkpoint_names=()

for fold_dir in "${fold_dirs[@]}"; do
    fold_base="$(basename "$fold_dir")"
    [[ "$fold_base" == fold_* ]] || exit_with_error "Invalid fold directory name: $fold_base"

    fold_id="${fold_base#fold_}"
    [[ -n "$fold_id" ]] || exit_with_error "Empty fold ID in directory name: $fold_base"
    fold_ids+=("$fold_id")

    mapfile -t checkpoint_files < <(find "$fold_dir" -mindepth 1 -maxdepth 1 -type f -name 'checkpoint_*.pth' | sort)
    [[ "${#checkpoint_files[@]}" -eq 1 ]] || exit_with_error "Expected exactly one checkpoint_* file in $fold_dir, found ${#checkpoint_files[@]}"

    checkpoint_base="$(basename "${checkpoint_files[0]}")"
    checkpoint_names+=("$checkpoint_base")
done

checkpoint_name="${checkpoint_names[0]}"
for chk in "${checkpoint_names[@]}"; do
    [[ "$chk" == "$checkpoint_name" ]] || exit_with_error "Not all folds have the same checkpoint filename: ${checkpoint_names[*]}"
done

#if [[ "${#fold_ids[@]}" -eq 1 ]]; then folds_arg="${fold_ids[0]}"; else folds_arg="($(IFS=,; echo "${fold_ids[*]}"))"; fi
if [[ "$use_gpu" == "1" ]]; then device="cuda"; else device="cpu"; fi

#make temporary directory and link input file to it
input_rt=$(basename "$input_file")
tempin=$(mktemp -d temp_in.XXX)
tempout=$(mktemp -d temp_out.XXX)
trap 'rm -r "$tempin" "$tempout"' EXIT
(cd $tempin && ln -sf "$input_file" "$input_rt" )

cmd=(
    nnUNetv2_predict
    -i "$tempin"
    -o "$tempout"
    -d "$dataset_id"
    -p "$plans_id"
    -tr "$trainer_id"
    -c "$configuration_id"
    -f "${fold_ids[@]}"
    -chk "$checkpoint_name"
    -npp "$nproc"
    -nps "$nproc"
    -device "$device"
)
#--disable_progress_bar

printf ' %q' "${cmd[@]}"; printf '\n'

"${cmd[@]}"


mapfile -t out_files < <( find "$tempout" -mindepth 1 -maxdepth 1 -type f \( -iname '*.nii' -o -iname \*.nii.gz \) )
[[ "${#out_files[@]}" -eq 1 ]] || exit_with_error "Expected exactly one output (.nii or .nii.gz) file, found ${#out_files[@]}"

#move output file
mv "${out_files[0]}" "$output_file"

#also move the remaining output files for provenance
mapfile -t out_files < <( find $tempout -maxdepth 1 -type f )
for file in "${out_files[@]}"; do b=$(basename "$file"); mv "$tempout/$b" "${output_file}-$b"; done
