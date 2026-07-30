#!/bin/bash

usage() {
    echo "Usage: $0 [-j <max_jobs>] [-q] [-s] [-T <timestep>] [-N <numsteps>] <advy_values> <h_values> <n_values> <p_values> [<method_values>]"
    echo "  -j max_jobs:     max parallel solver jobs (default: 1)"
    echo "  -q:              suppress solver output (still prints which run is starting)"
    echo "  -s:              skip runs whose .status file already exists in results_stability/"
    echo "  -T timestep:     TimeStep value (default: 0.01)"
    echo "  -N numsteps:     NumSteps value per run (default: 20)"
    echo "  advy_values:     comma-separated advy values        (e.g. \"0.5,1,2,4,8\")"
    echo "  h_values:        comma-separated AdaptBL_h_init     (e.g. \"0.2,0.1,0.05\")"
    echo "  n_values:        comma-separated NumRuns values     (e.g. \"1,16\")"
    echo "  p_values:        comma-separated NUMMODES values    (e.g. \"3,5,7\")"
    echo "  method_values:   comma-separated transfer methods   (default: \"Projection\")"
    echo
    echo "Runs that explode (NaN abort) do not stop the sweep; their solver exit code is"
    echo "recorded in the matching .status file alongside the .err file."
}

max_jobs=1
quiet=0
skip_existing=0
timestep=0.01
numsteps=20
while getopts "j:qsT:N:" opt; do
    case $opt in
        j) max_jobs="$OPTARG" ;;
        q) quiet=1 ;;
        s) skip_existing=1 ;;
        T) timestep="$OPTARG" ;;
        N) numsteps="$OPTARG" ;;
        *) usage; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

if [ "$#" -lt 4 ] || [ "$#" -gt 5 ]; then
    usage
    exit 1
fi

advy_values=$1
h_values=$2
n_values=$3
p_values=$4
method_values="${5:-Projection}"

if [ ! -f "ADR_stability_tmp.xml" ]; then
    echo "Error: ADR_stability_tmp.xml not found."
    exit 1
fi

IFS=',' read -ra advy_arr   <<< "$advy_values"
IFS=',' read -ra h_arr      <<< "$h_values"
IFS=',' read -ra n_arr      <<< "$n_values"
IFS=',' read -ra p_arr      <<< "$p_values"
IFS=',' read -ra method_arr <<< "$method_values"

resdir="results_stability"
mkdir -p "$resdir"

# Filename-safe form of a float: 0.05 -> 0dot05, -2.5 -> m2dot5
sanitise() {
    local v="${1//./dot}"
    echo "${v/#-/m}"
}

count=0
skipped=0
for advy in "${advy_arr[@]}"; do
    for h in "${h_arr[@]}"; do
        for n in "${n_arr[@]}"; do
            for p in "${p_arr[@]}"; do
                for method in "${method_arr[@]}"; do

                    advy_san=$(sanitise "$advy")
                    h_san=$(sanitise "$h")
                    stem="$resdir/ErrorFile_v_${advy_san}_h_${h_san}_n_${n}_p_${p}_m_${method}"
                    newfile="$resdir/ADR_v${advy_san}_h${h_san}_n${n}_p${p}_m${method}.xml"
                    log="$resdir/log_v${advy_san}_h${h_san}_n${n}_p${p}_m${method}.txt"

                    if [ "$skip_existing" -eq 1 ] && [ -f "${stem}.status" ]; then
                        (( skipped++ ))
                        continue
                    fi

                    cp ADR_stability_tmp.xml "$newfile"

                    sed -i "s/TRANSFER_METHOD/${method}/g" "$newfile"
                    if [[ "$method" == "ALE" ]]; then
                        sed -i 's|MOVEMENT_BLOCK|  <MOVEMENT>\n    <ZONES>\n      <CALLBACK ID="0" DOMAIN="D[0]" />\n    </ZONES>\n  </MOVEMENT>|' "$newfile"
                    else
                        sed -i '/MOVEMENT_BLOCK/d' "$newfile"
                    fi
                    sed -i "s/OUTPUT_METHOD/${method}/g"        "$newfile"
                    sed -i "s/ADVY_SAN/${advy_san}/g"           "$newfile"
                    sed -i "s/ADVY/${advy}/g"                   "$newfile"
                    sed -i "s/ADAPTBL_H_SAN/${h_san}/g"         "$newfile"
                    sed -i "s/ADAPTBL_H/${h}/g"                 "$newfile"
                    sed -i "s/NUMRUNS/${n}/g"                   "$newfile"
                    sed -i "s/NUMMODES_VAL/${p}/g"              "$newfile"
                    sed -i "s/TIMESTEP_VAL/${timestep}/g"       "$newfile"
                    sed -i "s/NUMSTEPS_VAL/${numsteps}/g"       "$newfile"

                    echo
                    echo "=== Starting: advy=${advy}, h_init=${h}, NumRuns=${n}, NUMMODES=${p}, method=${method} ==="
                    if [ "$quiet" -eq 0 ]; then
                        echo "  $NK1/ADRSolver-g squarecols.xml $newfile --force-output"
                    fi

                    (
                        # Never let an exploding run abort the sweep: capture the exit
                        # code and record it next to the .err file.
                        $NK1/ADRSolver-g squarecols.xml "$newfile" --force-output > "$log" 2>&1
                        rc=$?
                        echo "$rc" > "${stem}.status"
                        if [ "$quiet" -eq 1 ]; then
                            grep -i "error" "$log" | grep -iv "warning" || true
                        else
                            cat "$log"
                        fi
                        if [ "$rc" -ne 0 ]; then
                            echo "!!! EXPLODED (exit $rc): advy=${advy}, h_init=${h}, n=${n}, p=${p}, m=${method}"
                        fi
                    ) &

                    (( count++ ))
                    if (( count % max_jobs == 0 )); then
                        wait
                    fi

                done
            done
        done
    done
done

wait

echo
echo "=== Sweep summary ==="
ok=0
bad=0
for f in "$resdir"/*.status; do
    [ -e "$f" ] || continue
    if [ "$(cat "$f")" == "0" ]; then (( ok++ )); else (( bad++ )); fi
done
echo "Launched this invocation: $count"
if [ "$skip_existing" -eq 1 ]; then
    echo "Skipped (existing .status):  $skipped"
fi
echo "Total in $resdir/ - completed: $ok, exploded: $bad"
