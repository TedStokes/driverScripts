#!/bin/bash

if [ "$#" -ne 6 ]; then
    echo "Usage:"
    echo "bash fieldconvert_all.sh mesh.xml driver.xml chkprefix outprefix numruns numchksperrun"
    exit 1
fi

mesh_xml="$1"          # e.g. split1_v1.xml (v1 will be replaced automatically)
driver_xml="$2"
chk_prefix="$3"
out_prefix="$4"
numruns="$5"
numchks="$6"

start=0

for ((run=1; run<=numruns; run++)); do

    # Increment split1_vN.xml
    current_mesh=$(echo "$mesh_xml" | sed -E "s/v[0-9]+/v${run}/")

    if (( run == 1 )); then
        end=$(( start + numchks ))
    else
        end=$(( start + numchks + 1 ))
    fi

    echo "Run $run: converting $start → $end using $current_mesh"

    bash fieldconvert_chks.sh \
        "$current_mesh" \
        "$driver_xml" \
        "$chk_prefix" \
        "$out_prefix" \
        "$start" \
        "$end"

    start=$(( end + 1 ))

done