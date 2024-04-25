version 1.0

## Copyright Broad Institute, 2017
##
## This WDL workflow runs GATK4 Funcotator on a single tumor-normal pair or on a single tumor sample,
## and performs additional filtering and functional annotation tasks.
##
## Main requirements/expectations :
## - One analysis-ready BAM file (and its index) for each sample
##
## Description of inputs:
##
## ** Runtime **
## gatk_docker: docker image to use for GATK 4 Funcotator
## preemptible: how many preemptions to tolerate before switching to a non-preemptible machine (on Google)
## max_retries: how many times to retry failed tasks -- very important on the cloud when there are transient errors
## gatk_override: (optional) local file or Google bucket path to a GATK 4 java jar file to be used instead of the GATK 4 jar
##                in the docker image.  This must be supplied when running in an environment that does not support docker
##                (e.g. SGE cluster on a Broad on-prem VM)
##
## ** Workflow options **
## intervals optional): genomic intervals (will be used for scatter)
## funcotator_extra_args (optional): additional arguments for Funcotator calling and filtering
##
## ** Primary inputs **
## ref_fasta, ref_fai, ref_dict: reference genome, index, and dictionary
##
## Funcotator parameters (see Funcotator help for more details).
## funco_reference_version: "hg19" for hg19 or b37.  "hg38" for hg38.  Default: "hg19"
## funco_output_format: "MAF" to produce a MAF file, "VCF" to procude a VCF file.  Default: "MAF"
## funco_compress: (Only valid if funco_output_format == "VCF" )  If true, will compress the output of Funcotator.  If false, produces an uncompressed output file.  Default: false
## funco_use_gnomad_AF: If true, will include gnomAD allele frequency annotations in output by connecting to the internet to query gnomAD (this impacts performance).  If false, will not annotate with gnomAD.  Default: false
## funco_transcript_selection_mode: How to select transcripts in Funcotator.  ALL, CANONICAL, or BEST_EFFECT
## funco_transcript_selection_list: Transcripts (one GENCODE ID per line) to give priority during selection process.
## funco_data_sources_tar_gz:  Funcotator datasources tar gz file.  Bucket location is recommended when running on the cloud.
## funco_annotation_defaults:  Default values for annotations, when values are unspecified.  Specified as  <ANNOTATION>:<VALUE>.  For example:  "Center:Broad"
## funco_annotation_overrides:  Values for annotations, even when values are unspecified.  Specified as  <ANNOTATION>:<VALUE>.  For example:  "Center:Broad"
## funcotator_excluded_fields:  Annotations that should not appear in the output (VCF or MAF).  Specified as  <ANNOTATION>.  For example:  "ClinVar_ALLELEID"
## funcotator_extra_args: Any additional arguments to pass to Funcotator.  Default: ""
##
## Outputs :
## - One VCF/MAF file and its index with primary filtering applied; secondary filtering and functional annotation if requested; a bamout.bam
##   file of reassembled reads if requested
##
## Cromwell version support
## - Successfully tested on v34
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## pages at https://hub.docker.com/r/broadinstitute/* for detailed licensing information
## pertaining to the included programs.

struct Runtime {
    String gatk_docker
    File? gatk_override
    Int max_retries
    Int preemptible
    Int cpu
    Int machine_mem
    Int command_mem
    Int disk
    Int boot_disk_size
}

workflow funcotator_workflow {
    input {
      # Main inputs
      File ref_fasta
      File ref_fai
      File ref_dict
      File input_vcf
      File? input_vcf_idx

      String pair_id

      # Funcotator inputs
      String? case_id
      String? control_id
      String? sequencing_center
      String? sequence_source
      String? funco_reference_version
      String funco_output_format = "MAF"
      Boolean? funco_use_gnomad_AF
      File? funco_data_sources_tar_gz
      String? funco_transcript_selection_mode
      File? funco_transcript_selection_list
      Array[String]? funco_annotation_defaults
      Array[String]? funco_annotation_overrides
      Array[String]? funcotator_excluded_fields
      String? funcotator_extra_args

      String output_add = ".funcotated"
      Boolean compress = true

      # runtime
      String gatk_docker
      File? gatk_override
      Boolean filter_funcotations = true

      Int? preemptible
      Int? max_retries
      Int small_task_cpu = 2
      Int small_task_mem = 4
      Int small_task_disk = 100
      Int boot_disk_size = 12

      # Use as a last resort to increase the disk given to every task in case of ill behaving data
      Int? emergency_extra_disk

      # These are multipliers to multipler inputs by to make sure we have enough disk to accommodate for possible output sizes
      Float io_multiplier = 2.25
    }

    Int preemptible_or_default = select_first([preemptible, 2])
    Int max_retries_or_default = select_first([max_retries, 2])

    # Disk sizes used for dynamic sizing
    Int ref_size = ceil(size(ref_fasta, "GB") + size(ref_dict, "GB") + size(ref_fai, "GB"))
    Int input_vcf_size = ceil(size(input_vcf, "GB") + if defined(input_vcf_idx) then size(input_vcf_idx, "GB") else 0)

    # If no tar is provided, the task downloads one from broads ftp server
    Int funco_tar_size = if defined(funco_data_sources_tar_gz) then ceil(size(funco_data_sources_tar_gz, "GB") * 3) else 100
    Int gatk_override_size = if defined(gatk_override) then ceil(size(gatk_override, "GB")) else 0

    # This is added to every task as padding, should increase if systematically you need more disk for every call
    Int disk_pad = 10 + gatk_override_size + select_first([emergency_extra_disk,0])

    # Dynamic disk space
    Int disk_space = ceil(ref_size + input_vcf_size * io_multiplier + funco_tar_size + disk_pad)
    Int disk_space_or_default = select_first([small_task_disk, disk_space])

    # logic about output file names -- these are the names *without* .vcf extensions
    String output_name = sub(basename(input_vcf), "\\.vcf.*$", "") + output_add

    Runtime standard_runtime = {"gatk_docker": gatk_docker, "gatk_override": gatk_override,
            "max_retries": max_retries_or_default, "preemptible": preemptible_or_default, "cpu": small_task_cpu,
            "machine_mem": small_task_mem * 1000, "command_mem": small_task_mem * 1000 - 500,
            "disk": disk_space_or_default, "boot_disk_size": boot_disk_size}

    call Funcotate {
        input:
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            input_vcf = input_vcf,
            input_vcf_idx = input_vcf_idx,
            reference_version = select_first([funco_reference_version, "hg19"]),
            output_file_base_name = output_name,
            output_format = funco_output_format,
            compress = compress,
            use_gnomad = if defined(funco_use_gnomad_AF) then select_first([funco_use_gnomad_AF]) else false,
            data_sources_tar_gz = funco_data_sources_tar_gz,
            case_id = case_id,
            control_id = control_id,
            sequencing_center = sequencing_center,
            sequence_source = sequence_source,
            transcript_selection_mode = funco_transcript_selection_mode,
            transcript_selection_list = funco_transcript_selection_list,
            annotation_defaults = funco_annotation_defaults,
            annotation_overrides = funco_annotation_overrides,
            funcotator_excluded_fields = funcotator_excluded_fields,
            filter_funcotations = filter_funcotations,
            extra_args = funcotator_extra_args,
            runtime_params = standard_runtime
    }

    output {
        
        # Output files into Workspace
        File? funcotated_file = Funcotate.funcotated_output_file
        File? funcotated_file_index = Funcotate.funcotated_output_file_index
    }
}

task Funcotate {
     input {
       File ref_fasta
       File ref_fai
       File ref_dict
       File input_vcf
       File? input_vcf_idx
       String reference_version
       String output_file_base_name
       String output_format
       Boolean compress
       Boolean use_gnomad
       # This should be updated when a new version of the data sources is released
       # TODO: Make this dynamically chosen in the command.
       File? data_sources_tar_gz = "gs://broad-public-datasets/funcotator/funcotator_dataSources.v1.6.20190124s.tar.gz"
       String? control_id
       String? case_id
       String? sequencing_center
       String? sequence_source
       String? transcript_selection_mode
       File? transcript_selection_list
       Array[String]? annotation_defaults
       Array[String]? annotation_overrides
       Array[String]? funcotator_excluded_fields
       Boolean? filter_funcotations
       File? interval_list

       String? extra_args

       Runtime runtime_params
     }

     # ==============
     # Process input args:
     String output_maf = output_file_base_name + ".maf"
     String output_maf_index = output_maf + ".idx"
     String output_vcf = output_file_base_name + if compress then ".vcf.gz" else ".vcf"
     String output_vcf_idx = output_vcf +  if compress then ".tbi" else ".idx"
     String output_file = if output_format == "MAF" then output_maf else output_vcf
     String output_file_index = if output_format == "MAF" then output_maf_index else output_vcf_idx
     String transcript_selection_arg = if defined(transcript_selection_list) then " --transcript-list " else ""
     String annotation_def_arg = if defined(annotation_defaults) then " --annotation-default " else ""
     String annotation_over_arg = if defined(annotation_overrides) then " --annotation-override " else ""
     String filter_funcotations_args = if defined(filter_funcotations) && (filter_funcotations) then " --remove-filtered-variants " else ""
     String excluded_fields_args = if defined(funcotator_excluded_fields) then " --exclude-field " else ""
     String interval_list_arg = if defined(interval_list) then " -L " else ""
     String extra_args_arg = select_first([extra_args, ""])

     String dollar = "$"

     parameter_meta{
      ref_fasta: {localization_optional: true}
      ref_fai: {localization_optional: true}
      ref_dict: {localization_optional: true}
      input_vcf: {localization_optional: true}
      input_vcf_idx: {localization_optional: true}
     }

     command <<<
         set -e
         export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

         if [[ -f "~{input_vcf_idx}" ]]; then
             echo "VCF index exsits."
         else
             echo "Creating VCF index."
             gatk --java-options "-Xmx~{runtime_params.command_mem}m" IndexFeatureFile \
              -I ~{input_vcf}
         fi

         # Extract our data sources:
         echo "Extracting data sources zip file..."
         mkdir datasources_dir
         tar zxvf ~{data_sources_tar_gz} -C datasources_dir --strip-components 1
         DATA_SOURCES_FOLDER="$PWD/datasources_dir"

         # Handle gnomAD:
         if ~{use_gnomad} ; then
             echo "Enabling gnomAD..."
             for potential_gnomad_gz in gnomAD_exome.tar.gz gnomAD_genome.tar.gz ; do
                 if [[ -f ~{dollar}{DATA_SOURCES_FOLDER}/~{dollar}{potential_gnomad_gz} ]] ; then
                     cd ~{dollar}{DATA_SOURCES_FOLDER}
                     tar -zvxf ~{dollar}{potential_gnomad_gz}
                     cd -
                 else
                     echo "ERROR: Cannot find gnomAD folder: ~{dollar}{potential_gnomad_gz}" 1>&2
                     false
                 fi
             done
         fi

         # Run Funcotator:
         gatk --java-options "-Xmx~{runtime_params.command_mem}m" Funcotator \
             --data-sources-path $DATA_SOURCES_FOLDER \
             --ref-version ~{reference_version} \
             --output-file-format ~{output_format} \
             -R ~{ref_fasta} \
             -V ~{input_vcf} \
             -O ~{output_file} \
             ~{interval_list_arg} ~{default="" interval_list} \
             --annotation-default normal_barcode:~{default="Unknown" control_id} \
             --annotation-default tumor_barcode:~{default="Unknown" case_id} \
             --annotation-default Center:~{default="Unknown" sequencing_center} \
             --annotation-default source:~{default="Unknown" sequence_source} \
             ~{"--transcript-selection-mode " + transcript_selection_mode} \
             ~{transcript_selection_arg}~{default="" sep=" --transcript-list " transcript_selection_list} \
             ~{annotation_def_arg}~{default="" sep=" --annotation-default " annotation_defaults} \
             ~{annotation_over_arg}~{default="" sep=" --annotation-override " annotation_overrides} \
             ~{excluded_fields_args}~{default="" sep=" --exclude-field " funcotator_excluded_fields} \
             ~{filter_funcotations_args} \
             ~{extra_args_arg}
         # Make sure we have a placeholder index for MAF files so this workflow doesn't fail:
         if [[ "~{output_format}" == "MAF" ]] ; then
            touch ~{output_maf_index}
         fi
     >>>

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

     output {
         File funcotated_output_file = "~{output_file}"
         File funcotated_output_file_index = "~{output_file_index}"
     }
}