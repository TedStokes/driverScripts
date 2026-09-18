#!/bin/bash

usage() {
    echo "Usage: $0 [-j <max_jobs>] [-q] [-s] [-T <timestep>] [-N <numsteps>] [-P <x|y|xy>] [-V] [-o <dir>] <advy_values> <h_values> <n_values> <p_values> [<method_values>]"
    echo "  -j max_jobs:     max parallel solver jobs (default: 1)"
    echo "  -q:              suppress solver output (still prints which run is starting)"
    echo "  -s:              skip runs whose .status file already exists in the results dir"
    echo "  -T timestep:     TimeStep value (default: 0.01)"
    echo "  -N numsteps:     NumSteps value per run (default: 20)"
    echo "  -P dirs:         make the domain periodic in x, y or xy (default: Dirichlet"
    echo "                   exact-solution BCs on all sides). Uses squarecols_periodic.xml,"
    echo "                   whose boundary composites are ordered so that opposite edges"
    echo "                   pair up (see the top of this script). Filenames do not encode"
    echo "                   this either, so use a separate -o dir."
    echo "  -V:              also write a high-order vtu per timestep via a FieldConvert"
    echo "                   filter, into <dir>/vtus/sol_<run>_<step>_fc.vtu. Slow and"
    echo "                   bulky for long sweeps; meant for eyeballing single runs."
    echo "  -o dir:          results directory (default: results_stability). Filenames do"
    echo "                   not encode -T/-N, so give a separate dir when changing them"
    echo "                   or the previous sweep's results are overwritten."
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
periodic=""
fcfilter=0
resdir="results_stability"
while getopts "j:qsT:N:P:Vo:" opt; do
    case $opt in
        j) max_jobs="$OPTARG" ;;
        q) quiet=1 ;;
        s) skip_existing=1 ;;
        T) timestep="$OPTARG" ;;
        N) numsteps="$OPTARG" ;;
        P) periodic="$OPTARG" ;;
        V) fcfilter=1 ;;
        o) resdir="$OPTARG" ;;
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

case "$periodic" in
    ""|x|y|xy|yx) ;;
    *) echo "Error: -P must be one of x, y, xy (got '$periodic')."; exit 1 ;;
esac

if [ ! -f "ADR_stability_tmp.xml" ]; then
    echo "Error: ADR_stability_tmp.xml not found."
    exit 1
fi

# Nektar pairs periodic edges by their order within the boundary composite, and
# in squarecols.xml the top and bottom lists are in different x order - running
# periodic on it silently gives a wrong (much larger) error. The aligned mesh
# was made with the master-branch NekMesh ($NK1's build has no peralign module):
#   $MNK/NekMesh-g -m peralign:surf1=102:surf2=103:dir=y \
#                  -m peralign:surf1=104:surf2=105:dir=x \
#                  squarecols.xml squarecols_periodic.xml:xml:uncompress -f
mesh="squarecols.xml"
if [ -n "$periodic" ]; then
    mesh="squarecols_periodic.xml"
    if [ ! -f "$mesh" ]; then
        echo "Error: -P given but $mesh not found (regenerate it with NekMesh peralign;"
        echo "       see the comment in this script). squarecols.xml is NOT edge-aligned."
        exit 1
    fi
    echo "NOTE: periodic mode (-P $periodic) - using $mesh, not squarecols.xml,"
    echo "      because periodic edges must be paired in composite order."
fi

IFS=',' read -ra advy_arr   <<< "$advy_values"
IFS=',' read -ra h_arr      <<< "$h_values"
IFS=',' read -ra n_arr      <<< "$n_values"
IFS=',' read -ra p_arr      <<< "$p_values"
IFS=',' read -ra method_arr <<< "$method_values"

mkdir -p "$resdir"
if [ "$fcfilter" -eq 1 ]; then
    mkdir -p "$resdir/vtus"
fi

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

                    # Before any placeholder substitution, so RESDIR and the run
                    # placeholders inside the block get filled in too.
                    if [ "$fcfilter" -eq 1 ]; then
                        sed -i 's|FIELDCONVERT_BLOCK|    <FILTER TYPE="FieldConvert">\n      <PARAM NAME="OutputFile">RESDIR/vtus/sol_v_ADVY_SAN_h_ADAPTBL_H_SAN_n_NUMRUNS_p_NUMMODES_VAL_m_OUTPUT_METHOD.vtu:vtu:highorder</PARAM>\n      <PARAM NAME="OutputFrequency">1</PARAM>\n      <PARAM NAME="ResetCache">true</PARAM>\n    </FILTER>|' "$newfile"
                    else
                        sed -i '/FIELDCONVERT_BLOCK/d' "$newfile"
                    fi

                    # Must come first: the Error filter writes the .err file itself,
                    # so this is what actually honours -o. Piped delimiter because
                    # resdir may contain slashes.
                    sed -i "s|RESDIR|${resdir}|g"          "$newfile"
                    sed -i "s/TRANSFER_METHOD/${method}/g" "$newfile"
                    if [[ "$method" == "ALE" ]]; then
                        sed -i 's|MOVEMENT_BLOCK|  <MOVEMENT>\n    <ZONES>\n      <CALLBACK ID="0" DOMAIN="D[0]" />\n    </ZONES>\n  </MOVEMENT>|' "$newfile"
                    else
                        sed -i '/MOVEMENT_BLOCK/d' "$newfile"
                    fi
                    # Periodic BCs: swap each tagged Dirichlet line for a periodic
                    # condition pointing at the opposite region, then strip tags.
                    for dir in x y; do
                        if [[ "$periodic" == *"$dir"* ]]; then
                            DIR=${dir^^}
                            sed -i "s|^.*<!-- PERIODIC_${DIR}:\([0-9]*\) -->\$|        <P VAR=\"u\" VALUE=\"[\1]\" />|" "$newfile"
                        fi
                    done
                    sed -i 's| *<!-- PERIODIC_[XY]:[0-9]* -->||' "$newfile"
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
                    echo "=== Starting: advy=${advy}, h_init=${h}, NumRuns=${n}, NUMMODES=${p}, method=${method}${periodic:+, periodic=${periodic}} ==="
                    if [ "$quiet" -eq 0 ]; then
                        echo "  $NK1/ADRSolver-g $mesh $newfile --force-output"
                    fi

                    (
                        # Never let an exploding run abort the sweep: capture the exit
                        # code and record it next to the .err file.
                        $NK1/ADRSolver-g "$mesh" "$newfile" --force-output > "$log" 2>&1
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
