# Sourced by rundist to run in current session (no distribution)

snakemake -s ${SNAKEFILE} -j ${JOB_COUNT} \
    --ri -k \
    --jobname "{rulename}.{jobid}" \
    -w 30 "$@"
