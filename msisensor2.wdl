version 1.0

import "https://raw.githubusercontent.com/gatk-workflows/seq-format-conversion/master/cram-to-bam.wdl" as cram_to_bam

# Given a set of samples, combine segment files into a single file
# for tumor-only analysis
# The input tumor-only samples are expected as cram, which are then internally converted into bam and bai for msisensor2
workflow msisensor2_workflow {
    input {
        String pair_id
        File cram
        File crai

        # Inputs to go to cram_to_bam
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        
        Int memory=16
        Int disk_space = ceil((size(cram,"GB") + size(crai,"GB") + size(ref_fasta,"GB") + size(ref_fasta_index,"GB") + size(ref_dict,"GB")) * 4)
        Int num_threads = 20
        Int num_preempt = 3

    }

    call cram_to_bam.CramToBamFlow {
        input:
            ref_fasta=ref_fasta,
            ref_fasta_index=ref_fasta_index,
            ref_dict=ref_dict,
            input_cram=cram,
            sample_name=pair_id
    }

    call run_msisensor2 {
        input:
            sample_id=pair_id,
            bam=CramToBamFlow.outputBam,
            bai=CramToBamFlow.outputBai,
            memory=memory,
            disk_space=disk_space,
            num_threads=num_threads,
            num_preempt=num_preempt
    }

    output {
        File msisensor2_score = run_msisensor2.msisensor2_score
        File msisensor2_output = run_msisensor2.msisensor2_output
        File msisensor2_output_dis = run_msisensor2.msisensor2_output_dis
        File msisensor2_output_somatic = run_msisensor2.msisensor2_output_somatic
        File msisensor2_label = run_msisensor2.msisensor2_label
    }
}

task run_msisensor2 {

    input {
        String sample_id
        File bam
        File bai
        
        Int memory
        Int disk_space
        Int num_threads
        Int num_preempt
    }

    command <<<
        set -euo pipefail
    
        mv ~{bai} .
        mv ~{bam} .
        
        echo ~{bam} | rev | cut -d'/' -f1 | rev > new_bam_path.txt
        
        new_bam_path=$(cat new_bam_path.txt)
        msisensor2 msi -M /msisensor2/models_hg38 -t $new_bam_path -o ~{sample_id}.msisensor2.output
        head -2 ~{sample_id}.msisensor2.output | tail -1 | cut -f3 > ~{sample_id}.msisensor2.score

        # Creating label for exome sequenced sample
        python3 <<CODE
        with open("~{sample_id}.msisensor2.score",'r') as score_file:
            score = float(score_file.readlines()[0].split('\n')[0])
        
        if score >= 3.5:
            label="MSI"
        else:
            label="MSS"
        
        with open("~{sample_id}.msisensor2_label.txt", 'w') as label_file:
            label_file.write(label)
        CODE
    >>>

    output {
    	Float msisensor2_score=read_float("~{sample_id}.msisensor2.score")
        File msisensor2_output="~{sample_id}.msisensor2.output"
        File msisensor2_output_dis="~{sample_id}.msisensor2.output.dis"
        File msisensor2_output_somatic="~{sample_id}.msisensor2.output.somatic"
        File msisensor2_label="~{sample_id}.msisensor2_label.txt"
    }

    runtime {
        docker: "davidwu20/msisensor2:latest"
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
        preemptible: "~{num_preempt}"
    }
    
    meta {
        author: "Jivesh Singh"
    }
}