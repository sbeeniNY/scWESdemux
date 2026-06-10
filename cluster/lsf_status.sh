#!/usr/bin/env bash
set -u
# Map LSF job status to cluster-generic expected values: running / success / failed.

JOBID="${1:?Usage: lsf_status.sh JOBID}"
LOGFILE="logs/lsf/status_poll.log"
mkdir -p "$(dirname "${LOGFILE}")" 2>/dev/null || true

_log() {
  echo "[$(date '+%H:%M:%S')] job=${JOBID} $*" >> "${LOGFILE}"
}

status_from_bjobs() {
  bjobs -a -noheader -o 'stat' "${JOBID}" 2>/dev/null | awk 'NF { print $1; exit }'
}

status_from_bjobs_table() {
  bjobs -a "${JOBID}" 2>/dev/null | awk -v jobid="${JOBID}" 'NR > 1 && $1 == jobid { print $3; exit }'
}

status_from_history() {
  local history
  history="$(
    bhist -l "${JOBID}" 2>/dev/null ||
    bhist "${JOBID}" 2>/dev/null ||
    bacct -l "${JOBID}" 2>/dev/null ||
    true
  )"

  if grep -Eqi 'DONE|Done successfully|completed successfully' <<< "${history}"; then
    echo success
    return 0
  fi
  if grep -Eqi 'EXIT|Exited|exit code|TERM|Failed' <<< "${history}"; then
    echo failed
    return 0
  fi
  return 1
}

status_from_lsf_logs() {
  local err out
  err="$(find logs/lsf -type f -name "${JOBID}.e" -print -quit 2>/dev/null || true)"
  out="$(find logs/lsf -type f -name "${JOBID}.o" -print -quit 2>/dev/null || true)"

  if [[ -z "${err}${out}" ]]; then
    return 1
  fi

  if [[ -n "${out}" ]] && grep -Eqi 'Exited with exit code|TERM|Killed|Job was killed' "${out}"; then
    echo failed
    return 0
  fi

  if [[ -n "${err}" ]] && grep -Eqi 'Error in rule|WorkflowError|Traceback|segmentation fault|Permission denied|No such file|command not found' "${err}"; then
    echo failed
    return 0
  fi

  if [[ -n "${out}" ]] && grep -qi 'Successfully completed' "${out}"; then
    echo success
    return 0
  fi

  if [[ -n "${err}" ]] && grep -Eq 'Finished job [0-9]+\.|[0-9]+ of [0-9]+ steps \(100%\) done' "${err}"; then
    echo success
    return 0
  fi

  # LSF writes .o/.e files when a job finishes. If they exist but contain
  # no failure pattern, the job completed without errors.
  echo success
  return 0
}

# --- Main ---

STATUS="$(status_from_bjobs || true)"
SOURCE="bjobs-o"
if [[ -z "${STATUS}" ]]; then
  STATUS="$(status_from_bjobs_table || true)"
  SOURCE="bjobs-table"
fi

case "${STATUS}" in
  PEND|PSUSP|SSUSP|USUSP|WAIT)
    _log "source=${SOURCE} lsf=${STATUS} → running"
    echo running; exit 0 ;;
  RUN)
    _log "source=${SOURCE} lsf=${STATUS} → running"
    echo running; exit 0 ;;
  DONE)
    _log "source=${SOURCE} lsf=${STATUS} → success"
    echo success; exit 0 ;;
  EXIT|ZOMBI)
    _log "source=${SOURCE} lsf=${STATUS} → failed"
    echo failed;  exit 0 ;;
esac

HIST_RESULT="$(status_from_history || true)"
if [[ -n "${HIST_RESULT}" ]]; then
  _log "source=bhist/bacct → ${HIST_RESULT}"
  echo "${HIST_RESULT}"; exit 0
fi

LOG_RESULT="$(status_from_lsf_logs || true)"
if [[ -n "${LOG_RESULT}" ]]; then
  _log "source=lsf-logs → ${LOG_RESULT}"
  echo "${LOG_RESULT}"; exit 0
fi

_log "source=none (no bjobs/history/logs) → running"
echo running
