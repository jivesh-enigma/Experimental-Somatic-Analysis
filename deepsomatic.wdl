version 1.0

workflow DeepSomatic {
    input {
        String sample_id
        File ref_fasta
        File ref_fasta_index
        File cram
        File crai
        # Can be WGS,WES,PACBIO,ONT,FFPE_WGS,FFPE_WES,WGS_TUMOR_ONLY,PACBIO_TUMOR_ONLY,ONT_TUMOR_ONLY
        String DeepSomatic_Model_Type
        String permanent_bucket_path
        File interval_bed

        Int deepsomatic_ram_gb = 64
        Int deepsomatic_additional_disk_gb = 50
    }

    call deepsomatic_task {
        input:
            sample_id = sample_id,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            cram = cram,
            crai = crai,
            DeepSomatic_Model_Type = DeepSomatic_Model_Type,
            interval_bed = interval_bed,
            ram_gb = deepsomatic_ram_gb,
            additional_disk_gb = deepsomatic_additional_disk_gb
    }

    call deepsomatic_files_trasfer_permanent {
        input:
            sample_id = sample_id,
            deepsomatic_vcf = deepsomatic_task.deepsomatic_vcf,
            deepsomatic_vcf_index = deepsomatic_task.deepsomatic_vcf_index,
            deepsomatic_gvcf = deepsomatic_task.deepsomatic_gvcf,
            destination_bucket = permanent_bucket_path
    }

    output {
        String deepsomatic_vcf = deepsomatic_files_trasfer_permanent.deepsomatic_vcf_path
        String deepsomatic_vcf_index = deepsomatic_files_trasfer_permanent.deepsomatic_vcf_index_path
        String deepsomatic_gvcf = deepsomatic_files_trasfer_permanent.deepsomatic_gvcf_path
    }
}

task deepsomatic_task {
    
    input {
        String sample_id
        File ref_fasta
        File ref_fasta_index
        File cram
        File crai
        String DeepSomatic_Model_Type
        File interval_bed
        Int ram_gb
        Int additional_disk_gb
        Int input_size = ceil(size(cram,"GB") + size(ref_fasta,"GB")) * 4
        Int disk_gb = select_first([input_size, 100]) + additional_disk_gb
        Int preemptible = 1
        String DeepSomatic_docker = "google/deepsomatic:1.8.0"
    }

    command <<<
        mkdir intermediate_results_dir

        run_deepsomatic \
        --model_type=~{DeepSomatic_Model_Type} \
        --ref=~{ref_fasta} \
        --reads_tumor=~{cram} \
        --output_vcf=~{sample_id}.deepsomatic.vcf.gz \
        --output_gvcf=~{sample_id}.deepsomatic.g.vcf.gz \
        --regions=~{interval_bed} \
        --sample_name_tumor="~{sample_id}_DeepSomatic" \
        --num_shards=$(nproc) \
        --logging_dir=~{sample_id}.deepsomatic_logs \
        --intermediate_results_dir=intermediate_results_dir \
        --use_default_pon_filtering=true

        bcftools index -t ~{sample_id}.deepsomatic.vcf.gz
    >>>

    output {
        Array[File] deepsomatic_logs = glob("~{sample_id}.deepsomatic_logs/*")
        File deepsomatic_vcf = "~{sample_id}.deepsomatic.vcf.gz"
        File deepsomatic_vcf_index = "~{sample_id}.deepsomatic.vcf.gz.tbi"
        File deepsomatic_gvcf = "~{sample_id}.deepsomatic.g.vcf.gz"
    }

    runtime {
        docker: DeepSomatic_docker
        disks: "local-disk " + disk_gb + " HDD"
        memory: ram_gb + " GB"
        preemptible: preemptible
    }
}


task deepsomatic_files_trasfer_permanent {
    input {
        String sample_id
        String destination_bucket
        File deepsomatic_vcf
        File deepsomatic_vcf_index
        File deepsomatic_gvcf
    }

    command <<<
        gsutil -m cp ~{deepsomatic_vcf} ~{destination_bucket}/~{sample_id}/
        gsutil -m cp ~{deepsomatic_vcf_index} ~{destination_bucket}/~{sample_id}/
        gsutil -m cp ~{deepsomatic_gvcf} ~{destination_bucket}/~{sample_id}/
    >>>

    output {
        String deepsomatic_vcf_path = "~{destination_bucket}/~{sample_id}/~{sample_id}.deepsomatic.vcf.gz"
        String deepsomatic_vcf_index_path = "~{destination_bucket}/~{sample_id}/~{sample_id}.deepsomatic.vcf.gz.tbi"
        String deepsomatic_gvcf_path = "~{destination_bucket}/~{sample_id}/~{sample_id}.deepsomatic.g.vcf.gz"
    }

    runtime {
        docker: "google/cloud-sdk"
        memory: "4GB"
        disks: 'local-disk 50 HDD' 
        preemptible: "1"
    }
}