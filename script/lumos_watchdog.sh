#!/usr/bin/env bash
set -euo pipefail

# module load bioinfo-tools Nextflow samtools/1.20
# cd /proj/nobackup/sens2024549/wharf/shun-/shun--sens2024549/Lumos
# bash /proj/nobackup/sens2024549/Lumos/script/lumos_watchdog.sh

UPLOAD_ROOT="/proj/nobackup/sens2024549/wharf/shun-/shun--sens2024549/Lumos"
HVW_SCRIPT="/proj/nobackup/sens2024549/Lumos/script/run_lumos.sh"

CHECK_INTERVAL_SEC=300   # 5 min
STABLE_FOR_SEC=900       # 15 min
MAX_UPLOAD_SEC=$((24*3600))   # 24 hours hard timeout

STATE_DIR="${UPLOAD_ROOT}/.lumos_watchdog_state"

mkdir -p "$STATE_DIR"

PID_FILE="${STATE_DIR}/watchdog.pid"

if [[ -f "$PID_FILE" ]]; then
    old_pid="$(cat "$PID_FILE" 2>/dev/null || true)"

    if [[ -n "$old_pid" ]] && kill -0 "$old_pid" 2>/dev/null; then
        echo "Another watchdog is already running with PID ${old_pid}"
        exit 1
    fi

    rm -f "$PID_FILE"
fi

echo "$$" > "$PID_FILE"

WATCHDOG_LOG="${UPLOAD_ROOT}/watchdog.log"

rotate_log() {
    local logfile="$1"
    local max_lines=100

    [[ -f "$logfile" ]] || return 0

    local line_count
    line_count=$(wc -l < "$logfile")

    if (( line_count <= max_lines )); then
        return 0
    fi

    local tmpfile
    tmpfile="${logfile}.tmp"

    grep -Ev 'upload temp files detected|^(heartbeat\.)+$' "$logfile" \
		| sed -E 's/^(heartbeat\.)+//' \
		| tail -n "$max_lines" \
		> "$tmpfile"

    cat "$tmpfile" > "$logfile"
    rm -f "$tmpfile"
}

rotate_log "$WATCHDOG_LOG"

exec > >(tee -a "$WATCHDOG_LOG")
exec 2> >(tee -a "$WATCHDOG_LOG" >&2)

log() {
    printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

HPC_STATUS_LOG="${UPLOAD_ROOT}/HPC_status.log"

write_hpc_status() {
    {
        printf '# Time: [%s]\n' "$(date '+%F %T')"
        
        printf '\n# squeue\n'
        squeue -o "%.18i %.9P %.40j %.8u %.2t %.10M %.10l %.5C %.10m %.6D %R" || true

        printf '\n# projinfo 1 month\n' 
        projinfo || true
        
        printf '\n# projinfo 7 days\n'
        projinfo -s $(date -d '7 days ago' +%F) || true

        printf '\n# uquota\n'
        uquota || true
    } > "$HPC_STATUS_LOG"
}

cleanup() {
    trap - INT TERM EXIT

    log "Stopping watchdog and child monitoring jobs"

    rm -f "$PID_FILE"
    rm -f "$STATE_DIR"/*.uploading

    jobs -pr | xargs -r kill 2>/dev/null || true

    # Do not wait forever for background jobs.
    sleep 1
    jobs -pr | xargs -r kill -9 2>/dev/null || true

    exit 130
}

trap cleanup INT TERM
trap 'rm -f "$PID_FILE"' EXIT

bam_snapshot() {
    local dir="$1"

    find "$dir" -type f \
        \( -name '*.bam' -o -name '*.bai' -o -name '*.cram' -o -name '*.crai' \) \
        -printf '%p\t%s\n' 2>/dev/null \
    | sort
}

has_temp_files() {
    local dir="$1"

    find "$dir" -type f \
        \( -name '*.filepart' -o -name '*.part' -o -name '*.tmp' \) \
        -print -quit 2>/dev/null | grep -q .
}

sleep_interruptible() {
    sleep "$1" &
    wait $!
}

wait_until_stable() {
    local bam_dir="$1"
    local sample="$2"

    local start_time
    local now
    local total_elapsed

    local previous_snapshot=""
    local current_snapshot=""

    local stable_since=0
    local stable_elapsed=0

    start_time="$(date +%s)"

    log "${sample}: monitoring ${bam_dir}"

    while true; do
        now="$(date +%s)"
        total_elapsed=$(( now - start_time ))

        if (( total_elapsed > MAX_UPLOAD_SEC )); then
            log "${sample}: ERROR upload timeout"
            return 1
        fi

        if has_temp_files "$bam_dir"; then
            log "${sample}: upload temp files detected"

            stable_since=0
            previous_snapshot=""

            sleep "$CHECK_INTERVAL_SEC"
            continue
        fi

        current_snapshot="$(bam_snapshot "$bam_dir")"

        if [[ -z "$current_snapshot" ]]; then
            log "${sample}: no BAM-related files yet"

            stable_since=0
            previous_snapshot=""

            sleep "$CHECK_INTERVAL_SEC"
            continue
        fi

        if [[ "$current_snapshot" != "$previous_snapshot" ]]; then
            log "${sample}: BAM-related file list/size changed; waiting for stability"

            previous_snapshot="$current_snapshot"
            stable_since="$now"

            sleep "$CHECK_INTERVAL_SEC"
            continue
        fi

        stable_elapsed=$(( now - stable_since ))

        log "${sample}: BAM-related files stable for ${stable_elapsed}s"

        if (( stable_elapsed >= STABLE_FOR_SEC )); then
            log "${sample}: upload completed"
            return 0
        fi

        sleep_interruptible "$CHECK_INTERVAL_SEC"
    done
}

submit_workflow() {
    local sample="$1"
    local sample_dir="$2"
    local bam_dir="$3"

    local copied_script
    copied_script="${sample_dir}/$(basename "$HVW_SCRIPT")"

    log "${sample}: copying workflow script"
    cp -f "$HVW_SCRIPT" "$copied_script"
    # chmod +x "$copied_script"

    log "${sample}: submitting sbatch"

    sbatch \
        -A sens2024549 \
        -n 3 \
        -t 24:00:00 \
        -o "${sample_dir}/lumos-%j.out" \
        --wrap="cd '${sample_dir}' && bash '${copied_script}'"
}

log "Starting LUMOS watchdog"
log "UPLOAD_ROOT=${UPLOAD_ROOT}"

while true; do

    active_monitoring=0

    while read -r bam_dir; do

        sample_dir="$(dirname "$bam_dir")"
        sample="$(basename "$sample_dir")"

        find "$bam_dir" -type f -name '*.bam' -print -quit | grep -q . || continue

        done_file="${STATE_DIR}/${sample}.done"
        lock_file="${STATE_DIR}/${sample}.uploading"

        if [[ -f "$done_file" ]]; then
            continue
        fi

        if [[ -f "$lock_file" ]]; then
            active_monitoring=1
            continue
        fi

        touch "$lock_file"
        active_monitoring=1

        {
            trap 'rm -f "$lock_file"; exit 130' INT TERM

            if wait_until_stable "$bam_dir" "$sample"; then
                
                if [[ -f "$done_file" ]]; then
					log "${sample}: submission skipped; already marked done"
				else
					if submit_workflow "$sample" "$sample_dir" "$bam_dir"; then
						date '+%F %T' > "$done_file"
						log "${sample}: submitted successfully"
					else
						log "${sample}: ERROR sbatch submission failed"
						rm -f "$lock_file"
						exit 1
					fi
				fi
			else
				log "${sample}: monitoring failed"
			fi
        } &

    done < <(
        find "$UPLOAD_ROOT" \
            -mindepth 2 \
            -maxdepth 2 \
            -type d \
            -name bam
    )

    write_hpc_status

    if (( active_monitoring == 0 )); then
        # log "heartbeat: no active monitoring jobs"
        printf 'heartbeat.'
    fi

    sleep "$CHECK_INTERVAL_SEC"

done
