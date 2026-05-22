#!/bin/bash

usage() {
    echo "Usage: $0 [-j <max_jobs>] [-q] <io_checksteps> <n_values> <h_values> <p_values> [<method_values>]"
    echo "  -j max_jobs:     max parallel solver jobs (default: 1)"
    echo "  -q:              suppress solver output (still prints which run is starting)"
    echo "  -s:              skip runs whose .err file already exists in results/"
    echo "  n_values:        comma-separated NumRuns values    (e.g. \"1,2,4,8,16\")"
    echo "  h_values:        comma-separated AdaptBL_h values  (e.g. \"0.1,0.2,0.3\")"
    echo "  p_values:        comma-separated NUMMODES values   (e.g. \"3,5,7\")"
    echo "  method_values:   comma-separated transfer methods  (default: \"Projection\")"
}

max_jobs=1
quiet=0
skip_existing=0
while getopts "j:qs" opt; do
    case $opt in
        j) max_jobs="$OPTARG" ;;
        q) quiet=1 ;;
        s) skip_existing=1 ;;
        *) usage; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

if [ "$#" -lt 4 ] || [ "$#" -gt 5 ]; then
    usage
    exit 1
fi

io_checksteps=$1
n_values=$2
h_values=$3
p_values=$4
method_values="${5:-Projection}"

if ! [[ "$io_checksteps" =~ ^[0-9]+$ ]]; then
    echo "Error: io_checksteps must be a non-negative integer."
    exit 1
fi

if [ "$io_checksteps" -ne 0 ]; then
    echo "Error: io_checksteps output has been abandoned. Use VTU filter output instead."
    exit 1
fi

if [ ! -f "ADR_static_tmp.xml" ]; then
    echo "Error: ADR_static_tmp.xml not found."
    exit 1
fi

IFS=',' read -ra n_arr      <<< "$n_values"
IFS=',' read -ra h_arr      <<< "$h_values"
IFS=',' read -ra p_arr      <<< "$p_values"
IFS=',' read -ra method_arr <<< "$method_values"

mkdir -p results

# Suppress stdout/stderr of a command when -q is set
run_cmd() {
    if [ "$quiet" -eq 1 ]; then
        "$@" > /dev/null 2>&1
    else
        "$@"
    fi
}

count=0
skipped=0
for n in "${n_arr[@]}"; do
    for h in "${h_arr[@]}"; do
        for p in "${p_arr[@]}"; do
            for method in "${method_arr[@]}"; do

                h_san="${h//./dot}"
                method_san="${method}"
                newfile="results/ADR_n${n}_h${h_san}_p${p}_m${method_san}.xml"

                if [ "$skip_existing" -eq 1 ] && \
                   [ -f "results/ErrorFile_n_${n}_h_${h_san}_p_${p}_m_${method_san}.err" ]; then
                    (( skipped++ ))
                    continue
                fi

                cp ADR_static_tmp.xml "$newfile"

                sed -i "s/TRANSFER_METHOD/${method}/g" "$newfile"
                if [[ "$method" == "ALE" ]]; then
                    sed -i 's|MOVEMENT_BLOCK|  <MOVEMENT>\n    <ZONES>\n      <CALLBACK ID="0" DOMAIN="D[0]" />\n    </ZONES>\n  </MOVEMENT>|' "$newfile"
                else
                    sed -i '/MOVEMENT_BLOCK/d' "$newfile"
                fi
                sed -i "s/ADAPTBL_H_SAN/${h_san}/g"   "$newfile"
                sed -i "s/ADAPTBL_H/${h}/g"            "$newfile"
                sed -i "s/NUMRUNS/${n}/g"              "$newfile"
                sed -i "s/NUMMODES_VAL/${p}/g"         "$newfile"

                if [ "$io_checksteps" -eq 0 ]; then
                    sed -i "/IO_CHECKSTEPS/d" "$newfile"
                else
                    sed -i "s/IO_CHECKSTEPS/${io_checksteps}/g" "$newfile"
                fi

                echo
                echo "=== Starting: NumRuns=${n}, AdaptBL_h=${h}, NUMMODES=${p}, method=${method} ==="
                if [ "$quiet" -eq 0 ]; then
                    echo "  $NK1/ADRSolver-g square.xml $newfile --force-output"
                fi

                (
                    run_cmd $NK1/ADRSolver-g square.xml "$newfile" --force-output
                    if [ "$io_checksteps" -ne 0 ]; then
                        echo "bash fieldconvert_multi.sh square.xml $newfile square sol_n${n}_h${h_san}_p${p}_m${method_san} $n 1"
                        run_cmd bash fieldconvert_multi.sh square.xml "$newfile" square \
                            "sol_n${n}_h${h_san}_p${p}_m${method_san}" $n 1 $method_san
                    fi
                    # rm -f "$newfile"
                ) &

                (( count++ ))
                if (( count % max_jobs == 0 )); then
                    wait
                fi

            done
        done
    done
done

wait

if [ "$skip_existing" -eq 1 ]; then
    echo
    echo "Skipped $skipped run(s) with existing results."
fi
