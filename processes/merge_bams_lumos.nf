nextflow.enable.dsl = 2

process sortEachBam {
	tag "${params.sample}"
    cpus 1
    // memory '16 GB'
    time "30m"

    input:
    path bam

    output:
    path "${bam.simpleName}.primary.srt.bam"

    script:
    """
    set -euo pipefail

    samtools view -@ ${task.cpus} -h ${bam} \
      | awk '
          /^@/ { print; next }
          \$3 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)\$/ { print }
        ' \
      | samtools view -@ ${task.cpus} -b - \
      | samtools sort -@ ${task.cpus} -o "${bam.simpleName}.primary.srt.bam" -
    """
}

process mergeBamChunk {
	tag "${params.sample}"
    cpus 5
	memory "8 GB"
	time "6h"

    input:
    tuple val(idx), path(bams)

    output:
    tuple val(idx), path("chunk_${idx}.bam")

    script:
    def threads = Math.max(1, task.cpus - 2)
    """
    printf "%s\\n" ${bams} | tr ' ' '\\n' | sort > bam.list

    if [ \$(wc -l < bam.list) -eq 1 ]; then
        ln -s \$(cat bam.list) chunk_${idx}.bam
    else
        samtools merge -@ ${threads} -c -f -b bam.list chunk_${idx}.bam
    fi
    """
}

process mergeChunksFinal {
	tag "${params.sample}"
    cpus 7
    memory "8 GB"
    time "12h"
    publishDir "${params.outputDir}/bam", mode: "copy"

    input:
    path chunks

    output:
    tuple path("merged.bam"), path("merged.bam.bai")

    script:
    def threads = Math.max(1, task.cpus - 3)
    """
    find . -name 'chunk_*.bam' | sort > chunks.list

    if [ \$(wc -l < chunks.list) -eq 1 ]; then
        cp \$(cat chunks.list) merged.bam
        samtools index -@ ${threads} merged.bam
    else
        samtools merge -@ ${threads} -c -f \\
            -b chunks.list \\
            --write-index \\
            -o merged.bam##idx##merged.bam.bai
    fi
    """
}

process finalizeSingleChunk {
	tag "${params.sample}"
    cpus 2
    memory "8 GB"
    time "2h"
    publishDir "${params.outputDir}/bam", mode: "copy"

    input:
    path chunk

    output:
    tuple path("merged.bam"), path("merged.bam.bai")

    script:
    """
    cp ${chunk} merged.bam
    samtools index -@ ${task.cpus} merged.bam
    """
}

process indexSingleBam {
	tag "${params.sample}"
    publishDir "${params.outputDir}/bam", mode: "copy"

    cpus 4
    // memory "4 GB"
    time "1h"

    input:
    path bam

    output:
    tuple path("${bam.simpleName}.sorted.bam"), path("${bam.simpleName}.sorted.bam.bai")

    script:
    """
    samtools sort \
        -@ ${task.cpus} \
        --write-index \
        -o ${bam.simpleName}.sorted.bam##idx##${bam.simpleName}.sorted.bam.bai \
        ${bam}
    """
}

workflow ingress {
    take:
    bam_files

    main:
    //
    // collect all BAMs
    //
    bam_list_ch = bam_files.collect()

    //
    // SINGLE BAM
    //
    single_bam = bam_list_ch
        .filter { it.size() == 1 }
        .map { it[0] }

    single_with_index = single_bam
        .filter { file(it.toString() + ".bai").exists() }
        .map { bam ->
            tuple(
                bam,
                file(bam.toString() + ".bai")
            )
        }

    single_without_index = single_bam
        .filter { !file(it.toString() + ".bai").exists() }

    indexSingleBam(single_without_index)

    //
    // MULTI BAM
    //
    multi_bam = bam_list_ch
        .filter { it.size() > 1 }
        .flatMap { it }

    sorted = sortEachBam(multi_bam)

    chunked = sorted
        .collect()
        .flatMap { xs ->
            xs.collate(50)
              .withIndex()
              .collect { chunk, idx ->
                    tuple(idx, chunk)
              }
        }

    mergeBamChunk(chunked)

    all_chunks = mergeBamChunk.out
        .map { idx, bam -> bam }
        .collect()

    single_chunk = all_chunks
		.filter { it.size() == 1 }
		.map { it[0] }

	multi_chunks = all_chunks
		.filter { it.size() > 1 }

	finalizeSingleChunk(single_chunk)
	mergeChunksFinal(multi_chunks)

    //
    // unify outputs
    //
    final_pairs = single_with_index
        .mix(indexSingleBam.out)
        .mix(finalizeSingleChunk.out)
        .mix(mergeChunksFinal.out)

    bam_out = final_pairs.map { bam, bai -> bam }
    bai_out = final_pairs.map { bam, bai -> bai }

    emit:
    bam = bam_out
    bai = bai_out
}
