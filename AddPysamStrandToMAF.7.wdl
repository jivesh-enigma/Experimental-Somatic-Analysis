version 1.0

workflow AddPysamStrandToMAF {
    input {
        File tumor_bam
        File tumor_bai
        File maf
    }

    call AddStrandToMAF_Task {
        input:
            tumor_bam = tumor_bam,
            tumor_bai = tumor_bai,
            maf = maf
    }

    output {
        File strand_added_maf = AddStrandToMAF_Task.strand_added_maf
    }
}

task AddStrandToMAF_Task {
    input {
        File tumor_bam
        File tumor_bai
        File maf

        String outbase = basename(maf, ".maf")

        # Docker Configs
        String docker = "hoyinchudfci/pysam_add_strand_to_maf:1.0"
        Int diskSpaceGB = 256
        Int memoryGb = 4
        Int cpu = 1
        Int bootSpaceGB = 256
    }

    command {
        pwd
        ls -al
        
        cd ..
        
        echo "cd"
        pwd
        ls -al
        
        mv ~{tumor_bam} ./tumor.bam
        mv ~{tumor_bai} ./tumor.bai
        
        echo "moved"
        pwd
        ls -al

        echo "finding python file"
        echo "up"
        ls ..
        echo "UP"
        ls ../..

        python ./pysam_strand_check.py \
        -bam ./tumor.bam \
        -maf ~{maf} \
        -out ~{outbase}.strand_added.maf

        mv ~{outbase}.strand_added.maf /cromwell_root/
        
        ls -al
    }

    output {
        File strand_added_maf = "${outbase}.strand_added.maf"
    }

    runtime {
        docker: docker
        memory: "${memoryGb} GB"
        cpu: "${cpu}"
        disks: "local-disk ${diskSpaceGB} HDD"
        bootDiskSizeGb: bootSpaceGB
    }
}
