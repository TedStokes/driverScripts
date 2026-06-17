#!/bin/bash

if [ "$#" -ne 7 ]; then
    echo "Usage:"
    echo "bash fieldconvert_multi.sh mesh.xml driver.xml chkprefix outprefix numruns numsteps {ALE|Projection}"
    exit 1
fi

mesh_xml="$1"          # e.g. split1_v1.xml (v1 will be replaced automatically)
driver_xml="$2"
chk_prefix="$3"
out_prefix="$4"
numruns="$5"
numchks="$6"           # NumSteps (= IO_CheckSteps checkpoints per DoSolve)
transfer="$7"          # ALE or Projection

# Checkpoint counts per mesh version (N = numchks = NumSteps, assuming IO_CheckSteps=1).
# The initial DoSolve outputs N+1 checkpoints (step 0 initial condition + N steps).
#
# Projection:
#   run=1       : N+1 chks  (step-0 + N steps; pre-transfer goes into run=2's range)
#   run=2..last : N+2 chks  (pre-transfer 1 + post-proj 1 + DoSolve N)
#
# ALE:
#   run=1       : N+2 chks  (step-0 + N steps + pre-transfer 1)
#   run=2..n-1  : N+1 chks  (post-ALE 1 + N-1 remaining + pre-transfer 1)
#   run=n (last): N   chks  (post-ALE 1 + N-1 remaining, no next pre-transfer)

start=0

# For ALE, both the mesh XML and driver XML contain a <MOVEMENT> block with a
# CALLBACK zone.  FieldConvert never registers these callbacks (only ADRSolver
# does), so OutputVtk::v_OutputFromExpERN crashes when it calls PerformMovement
# on an empty std::function.  Strip <MOVEMENT> from both XMLs; the geometry is
# already baked into each .chk file at checkpoint time.
fc_driver_xml="$driver_xml"
tmp_driver_xml=""
strip_movement() { sed '/<MOVEMENT>/,/<\/MOVEMENT>/d' "$1"; }
if [[ "$transfer" == "ALE" ]]; then
    tmp_driver_xml=$(mktemp --suffix=.xml)
    strip_movement "$driver_xml" > "$tmp_driver_xml"
    fc_driver_xml="$tmp_driver_xml"
fi

for ((run=1; run<=numruns; run++)); do

    # Increment split1_vN.xml
    current_mesh=$(echo "$mesh_xml" | sed -E "s/v[0-9]+/v${run}/")

    if [[ "$transfer" == "ALE" ]]; then
        if (( run == 1 )); then
            end=$(( start + numchks + 1 ))    # initial DoSolve (N+1 chks incl. step-0) + pre-transfer
        elif (( run == numruns )); then
            end=$(( start + numchks - 1 ))    # post-ALE + N-1 remaining, no next pre-transfer
        else
            end=$(( start + numchks ))        # post-ALE + N-1 remaining + pre-transfer
        fi
    else  # Projection (unchanged from original)
        if (( run == 1 )); then
            end=$(( start + numchks ))
        else
            end=$(( start + numchks + 1 ))
        fi
    fi

    echo "Run $run: converting $start → $end using $current_mesh"

    fc_mesh_xml="$current_mesh"
    tmp_mesh_xml=""
    if [[ "$transfer" == "ALE" ]]; then
        tmp_mesh_xml=$(mktemp --suffix=.xml)
        strip_movement "$current_mesh" > "$tmp_mesh_xml"
        fc_mesh_xml="$tmp_mesh_xml"
    fi

    # echo "bash fieldconvert_chks.sh $fc_mesh_xml $fc_driver_xml $chk_prefix $out_prefix $start $end"
    bash fieldconvert_chks.sh \
        "$fc_mesh_xml" \
        "$fc_driver_xml" \
        "$chk_prefix" \
        "$out_prefix" \
        "$start" \
        "$end"

    [[ -n "$tmp_mesh_xml" ]] && rm -f "$tmp_mesh_xml"
    start=$(( end + 1 ))

done

[[ -n "$tmp_driver_xml" ]] && rm -f "$tmp_driver_xml"
