#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
converter="$script_dir/convert_csv"

if [[ ! -x $converter ]]; then
    echo "Converter executable not found: $converter" >&2
    echo "Build it with: make -C '$script_dir' convert_csv" >&2
    exit 1
fi

shopt -s nullglob
files=(field_config_SU*.csv)
matched=0

for file in "${files[@]}"; do
    if [[ ! $file =~ ^field_config_SU([0-9]+)_([0-9]+([.][0-9]+)?)_([0-9]+)x([0-9]+)x([0-9]+)(_[0-9]+)?\.csv$ ]]; then
        continue
    fi

    ncol=${BASH_REMATCH[1]}
    beta=${BASH_REMATCH[2]}
    ls_first=${BASH_REMATCH[4]}
    ls_second=${BASH_REMATCH[5]}
    lt=${BASH_REMATCH[6]}
    run=${BASH_REMATCH[7]:-}

    # The supported filename form uses the same LS in both spatial dimensions.
    [[ $ls_first == "$ls_second" ]] || continue

    args=(--ncol "$ncol" --beta "$beta" --ls "$ls_first" --lt "$lt")
    if [[ -n $run ]]; then
        args+=(--run "$run")
    fi

    printf 'Converting %s\n' "$file"
    "$converter" "${args[@]}"
    matched=$((matched + 1))
done

if ((matched == 0)); then
    echo "No matching configuration files found." >&2
    exit 1
fi
