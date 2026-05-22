#!/bin/bash
MAX_JOBS=11
count=0

for s in $(seq $5 $6); do
    # echo "$NK1/FieldConvert-g $1 $2 ${3}_${s}.chk ${4}_${s}.vtu:vtu:highorder -n 5 -f"
    $NK1/FieldConvert-g "$1" "$2" "${3}_${s}.chk" "${4}_${s}.vtu:vtu:highorder" -n 5 -f &
    ((count++))
    if (( count % MAX_JOBS == 0 )); then
        wait
    fi
done

wait
