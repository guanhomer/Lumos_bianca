# module load bioinfo-tools Nextflow samtools/1.20
# cd /proj/nobackup/sens2024549/Lumos
# sbatch -A sens2024549 -n 3 -t 24:00:00 --wrap="bash run_lumos.sh" -o /proj/nobackup/sens2024549/wharf/shun-/shun--sens2024549/slurm-%j.out

# Nextflow settings
#ls /sw/bioinfo/Nextflow/latest/rackham/nxf_home/framework/
export NXF_VER=25.10.0
export NXF_HOME=/castor/project/proj_nobackup/tools/nextflow/.nextflow
export NXF_SINGULARITY_CACHEDIR=/proj/nobackup/sens2024549/Lumos/KolmogorovLab/singularity
export NXF_OFFLINE=true
# export NXF_OPTS='-Xms1g -Xmx4g'

export TMPDIR="${SLURM_TMPDIR:-/scratch/$USER/tmp}"
mkdir -p "$TMPDIR"
export NXF_TEMP="$TMPDIR"
export NXF_OPTS="-Xms1g -Xmx4g -Djava.io.tmpdir=$TMPDIR"

# Workflow paths
WORKFLOW_DIR=/proj/nobackup/sens2024549/Lumos/KolmogorovLab/Lumos
CUSTOM_CONFIG=$WORKFLOW_DIR/uppmax_ryan.config

# update status
SAMPLE="$(basename "$PWD")"
STATE_DIR="/proj/nobackup/sens2024549/wharf/shun-/shun--sens2024549/Lumos/.lumos_watchdog_state"
PROCESSING_FILE="${STATE_DIR}/${SAMPLE}.processing"
mkdir -p "$STATE_DIR"
touch "$PROCESSING_FILE"
cleanup() {
    rm -f "$PROCESSING_FILE"
}
trap cleanup EXIT INT TERM

# Run workflow
cmd="nextflow run $WORKFLOW_DIR/tumorOnlyONT.nf"
#cmd+=" -profile test"
cmd+=" -c $CUSTOM_CONFIG"
cmd+=" --project sens2024549"
cmd+=" --reference \"/proj/sens2024549/reference/GRCh38.p14.genome.fa\""
cmd+=" --clair3_model \"${WORKFLOW_DIR}/clair3_models_pytorch/r1041_e82_400bps_sup_v520\""
cmd+=" --sv_pon \"${WORKFLOW_DIR}/annot/PoN_1000G_hg38_extended.tsv.gz\""
cmd+=" --vntr \"${WORKFLOW_DIR}/annot/human_GRCh38_no_alt_analysis_set.trf.bed\""
cmd+=" --cpgs \"${WORKFLOW_DIR}/annot/hg38_cpg_cleaned.bed\""
cmd+=" --cosmic \"${WORKFLOW_DIR}/annot/cosmic_genes.tsv\""
#cmd+=" --reads_tumor \"testdata/*.bam\""
cmd+=" --mode all"
#cmd+=" --no_dmr"
cmd+=" --aligned_input"
cmd+=" --aligned_tumor_dir bam"
cmd+=" --sample ${SAMPLE}"
cmd+=" -resume"

# Print for debug
echo "Running command:"
echo "$cmd"

# Execute
eval "$cmd"

