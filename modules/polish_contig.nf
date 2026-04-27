process polish_contig {
    label 'polish'
    
    tag "${sample_id}"

    publishDir "${params.outdir}/assembled_genomes", mode: 'copy', saveAs: { filename ->
        filename.endsWith("_polished.fasta") ? filename : null
    }

    input:
    tuple val(sample_id), path(assembly_fasta), path(r1_trimmed), path(r2_trimmed)

    output:
    path "${sample_id}_polished.fasta", emit: polished_assembly
    path "${sample_id}_polishing_report.txt", emit: polishing_report
    path "${sample_id}_polishing_bam.bam", emit: bam
    path "${sample_id}_polishing_bam.bam.bai", emit: bai

    script:
    """
    ############################
    # VARIABLES
    ############################
    ASSEMBLY="${assembly_fasta}"
    R1="${r1_trimmed}"
    R2="${r2_trimmed}"
    
    OUTDIR="pilon_workdir"
    mkdir -p "\$OUTDIR"

    ############################
    # STEP 1 — INDEX ASSEMBLY
    ############################
    bwa-mem2 index "\$ASSEMBLY"
    samtools faidx "\$ASSEMBLY"

    ############################
    # STEP 2 — READ MAPPING
    ############################
    BAM="\$OUTDIR/aligned.bam"

    bwa-mem2 mem -t ${params.cpus} -K ${params.pol_bwa_K} -Y "\$ASSEMBLY" "\$R1" "\$R2" -B ${params.pol_bwa_B} -L ${params.pol_bwa_L} | \\
        samclip --ref "\$ASSEMBLY.fai" | \\
        samtools sort -@ ${params.cpus} -o "\$BAM"

    samtools index "\$BAM"

    ############################
    # STEP 3 — PILON POLISHING
    ############################
    pilon \\
        --genome "\$ASSEMBLY" \\
        -Xmx${params.ram}G \\
        --frags "\$BAM" \\
        --output pilon \\
        --threads ${params.cpus} \\
        --fix all \\
        --changes \\
        --mindepth ${params.pol_mindepth} \\
        --minmq ${params.pol_minmq} \\
        --minqual ${params.pol_minqual} \\
        --mingap ${params.pol_mingap}

    ############################
    # STEP 4 — FINAL ASSEMBLY
    ############################
    mv pilon.fasta "${sample_id}_polished.fasta"
    mv pilon.changes pilon.changes.txt

    ############################
    # STEP 4B — EXTRACT ORIGINAL STATS
    ############################
    grep "^>" "\$ASSEMBLY" | awk '{
        header = \$1
        sub(/^>/, "", header)

        ori_cov = "NA"
        pos = index(header, "_cov_")
        if (pos > 0) {
            cov_str = substr(header, pos + 5)
            n = split(cov_str, parts, "_")
            ori_cov = parts[1]
        }

        print header "\\t" ori_cov
    }' > original_stats.txt

    ############################
    # STEP 5 — DEPTH PER CONTIG
    ############################
    samtools depth -a "\$BAM" | \\
        awk '{
            depth[\$1] += \$3
            count[\$1]++
        }
        END {
            for (c in depth)
                printf "%s\\t%.2f\\n", c, depth[c]/count[c]
        }' > per_contig_depth.txt

    ############################
    # STEP 6 — PARSE PILON CHANGES
    ############################
    TOTAL_CONTIGS=\$(grep -c "^>" "\$ASSEMBLY")

    awk -v total_contigs="\$TOTAL_CONTIGS" '{
        split(\$1, parts, ":")
        contig = parts[1]
        orig   = \$3
        newb   = \$4

        # Count the full bases (ATCGN)
        if (orig ~ /[Nn]/ && newb != ".") {
            gap_fills[contig]++
            gap_bases_filled[contig] += length(newb)
            gap_bases_original[contig] += length(orig)

        # Local rearrangement: both orig and new are long sequences (>10 bp)
        } else if (length(orig) > 10 && length(newb) > 10 && orig != newb) {
            local_fix[contig]++

        # Insertion: new bases added (orig is ".")
        } else if (orig == "." || length(newb) > length(orig)) {
            ins[contig] += length(newb)

        # Deletion: bases removed (new is ".")
        } else if (newb == "." || length(orig) > length(newb)) {
            dels[contig] += length(orig)

        # Substitution: single base change
        } else {
            subs[contig]++
        }

        total[contig]++
    }
    END {
        for (c in total) {
            s  = (subs[c]       ? subs[c]       : 0)
            i  = (ins[c]        ? ins[c]        : 0)
            d  = (dels[c]       ? dels[c]       : 0)
            gf = (gap_fills[c]  ? gap_fills[c]  : 0)
            gbo = (gap_bases_original[c] ? gap_bases_original[c] : 0)
            gbf = (gap_bases_filled[c]   ? gap_bases_filled[c]   : 0)
            lf = (local_fix[c]  ? local_fix[c]  : 0)
            printf "%s\\t%d\\t%d\\t%d\\t%d\\t%d\\t%d\\t%d\\t%d\\n", \
                c, total[c], s, i, d, gf, gbo, gbf, lf
        }
        printf "TOTAL_CONTIGS\\t%d\\t0\\t0\\t0\\t0\\t0\\t0\\t0\\n", total_contigs
    }' pilon.changes.txt > pilon_status.txt

    ############################
    # STEP 7 — ANNOTATE HEADERS
    ############################
    awk '
        BEGIN {
            while ((getline line < "per_contig_depth.txt") > 0) {
                split(line, a, "\\t")
                depth[a[1]] = a[2]
            }
            while ((getline line < "pilon_status.txt") > 0) {
                split(line, a, "\\t")
                total[a[1]] = a[2]
                subs[a[1]]  = a[3]
                ins[a[1]]   = a[4]
                dels[a[1]]  = a[5]
            }
            while ((getline line < "original_stats.txt") > 0) {
                split(line, a, "\\t")
                ori_cov[a[1]] = a[2]
            }
        }
        /^>/ {
            contig = substr(\$0, 2)
            display = contig
            if (substr(display, length(display)-5, 6) == "_pilon")
                display = substr(display, 1, length(display)-6)
            lookup = display
            d  = (lookup in depth) ? depth[lookup] : "0.00"
            t  = (lookup in total) ? total[lookup] : "0"
            s  = (lookup in subs)  ? subs[lookup]  : "0"
            i  = (lookup in ins)   ? ins[lookup]   : "0"
            dl = (lookup in dels)  ? dels[lookup]  : "0"
            oc = (lookup in ori_cov) ? ori_cov[lookup] : "NA"
            print ">" display "_ori_cov_" oc "_polish_mean_depth_" d "X_polished=" t "_sub=" s "_ins=" i "_del=" dl
            next
        }
        { print }
    ' "${sample_id}_polished.fasta" > tmp.fasta

    mv tmp.fasta "${sample_id}_polished.fasta"

    ############################
    # STEP 7C — SORT BY LENGTH & RENAME
    ############################
    # Calculate lengths and sort contigs by size (descending)
    awk '
        /^>/ {
            if (seq) {
                print header
                print seq
            }
            header = \$0
            seq = ""
            next
        }
        { seq = seq \$0 }
        END {
            if (seq) {
                print header
                print seq
            }
        }
    ' "${sample_id}_polished.fasta" | \\
    awk '
        /^>/ {
            if (name) {
                print name "\\t" length(sequence) "\\t" sequence
            }
            name = \$0
            sequence = ""
            next
        }
        { sequence = sequence \$0 }
        END {
            if (name) {
                print name "\\t" length(sequence) "\\t" sequence
            }
        }
    ' | \\
    sort -t\$'\\t' -k2,2nr | \\
    awk -F'\\t' '
        {
            rank++
            # Extract fields
            header = \$1
            new_length = \$2
            sequence = \$3
            
            # Remove leading ">"
            sub(/^>/, "", header)
            
            # Extract the original node name
            if (match(header, /^[^_]+_[0-9]+/)) {
                ori_id = tolower(substr(header, RSTART, RLENGTH))
            } else if (match(header, /^[^_]+/)) {
                ori_id = tolower(substr(header, RSTART, RLENGTH))
            } else {
                ori_id = tolower(header)
            }
            
            # Extract original coverage
            ori_cov = "NA"
            if (match(header, /_ori_cov_[^_]+/)) {
                cov_field = substr(header, RSTART)
                if (match(cov_field, /[0-9.]+/)) {
                    ori_cov = substr(cov_field, RSTART, RLENGTH)
                }
            }
            
            # Extract polishing annotations (from _polish_mean_depth_ onwards)
            polish_annotations = ""
            if (match(header, /_polish_mean_depth_.*/)) {
                polish_annotations = substr(header, RSTART)
            }
            
            # Create new header with new length and original coverage
            new_header = sprintf(">polished_node_%d_ori_%s_len_%d_cov_%s%s", 
                rank, ori_id, new_length, ori_cov, polish_annotations)
            print new_header
            print sequence
        }
    ' > "${sample_id}_polished_sorted.fasta"

    mv "${sample_id}_polished_sorted.fasta" "${sample_id}_polished.fasta"

    ############################
    # STEP 7B — COMPARE POLISHED
    ############################
    python3 ${projectDir}/bin/compare_polished.py "\$ASSEMBLY" "${sample_id}_polished.fasta"

    ############################
    # STEP 8 — EXPORT BAM
    ############################
    cp "\$BAM" "${sample_id}_polishing_bam.bam"
    cp "\$BAM.bai" "${sample_id}_polishing_bam.bam.bai"

############################
    # STEP 9 — REPORT
    ############################
    {
        echo "=========================================="
        echo "POLISHING REPORT"
        echo "=========================================="
        echo "Sample: ${sample_id}"
        echo "Date: \$(date)"
        echo ""

        echo "Assembly stats:"
        seqkit stats -T "${sample_id}_polished.fasta" | tail -n 1 | \
            awk '{printf "  Contigs: %s\\n  Total length: %s bp\\n", \$4, \$5}'

        echo ""
        echo "Mapping stats:"
        samtools flagstat "\$BAM" | head -5

        echo ""
        echo "Coverage:"
        samtools depth "\$BAM" | \
            awk '{sum+=\$3; count++} END {printf "  Mean depth: %.2f X\\n", sum/count}'

        echo ""
        echo "Correction summary:"
        awk '{
            if (\$1 == "TOTAL_CONTIGS") {
                all = \$2
            } else {
                events += \$2; subs += \$3; ins += \$4; dels += \$5
                total_bases += \$3 + \$4 + \$5
                if (\$2 > 0) corrected++
            }
        }
        END {
            printf "  Polished contigs: %d / %d\\n", corrected, all
            printf "  Total events: %d\\n", events
            printf "  Total bases changed: %d\\n", total_bases
            printf "  Substituted bases: %d\\n", subs
            printf "  Inserted bases: %d\\n", ins
            printf "  Deleted bases: %d\\n", dels
        }' pilon_status.txt

        echo ""
        echo "Gap filling summary:"
        awk '{
            if (\$1 != "TOTAL_CONTIGS" && \$6 > 0) {
                gaps_filled += \$6
                original_gap_bases += \$7
                bases_filled += \$8
                contigs_with_gaps++
            }
        }
        END {
            printf "  Contigs with gaps filled: %d\\n", contigs_with_gaps
            printf "  Total gaps filled: %d\\n", gaps_filled
            printf "  Total bases filled: %d/%d\\n", original_gap_bases, bases_filled
        }' pilon_status.txt

        echo ""
        echo "Local misassembly rearrangements:"
        awk '{
            if (\$1 != "TOTAL_CONTIGS" && \$9 > 0) {
                local_total += \$9
                contigs_rearranged++
            }
        }
        END {
            printf "  Contigs rearranged: %d\\n", contigs_rearranged
            printf "  Total local fixes: %d\\n", local_total
        }' pilon_status.txt

        echo ""
        echo "✅ Pilon polishing completed"
        echo "=========================================="
    } > "${sample_id}_polishing_report.txt"
    """
}
