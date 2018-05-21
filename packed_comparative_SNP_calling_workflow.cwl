{
    "cwlVersion": "v1.0", 
    "$graph": [
        {
            "class": "CommandLineTool", 
            "label": "index a compressed *.vcf.gz file", 
            "baseCommand": [
                "bcftools", 
                "index"
            ], 
            "requirements": [
                {
                    "class": "DockerRequirement", 
                    "dockerPull": "quay.io/biocontainers/bcftools:1.8--1"
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "class": "InitialWorkDirRequirement", 
                    "listing": [
                        "$(inputs.input_vcf)"
                    ]
                }
            ], 
            "inputs": [
                {
                    "id": "#bcftools_index.cwl/input_vcf", 
                    "type": "File", 
                    "inputBinding": {
                        "position": 1
                    }
                }
            ], 
            "outputs": [
                {
                    "id": "#bcftools_index.cwl/indexed_vcf", 
                    "type": "File", 
                    "secondaryFiles": [
                        ".csi"
                    ], 
                    "outputBinding": {
                        "glob": "$(inputs.input_vcf.basename)"
                    }
                }
            ], 
            "id": "#bcftools_index.cwl"
        }, 
        {
            "class": "CommandLineTool", 
            "label": "remove indels from vcf", 
            "baseCommand": [
                "bcftools", 
                "isec"
            ], 
            "requirements": [
                {
                    "class": "DockerRequirement", 
                    "dockerPull": "quay.io/biocontainers/bcftools:1.8--1"
                }, 
                {
                    "class": "StepInputExpressionRequirement"
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "class": "InitialWorkDirRequirement", 
                    "listing": [
                        "$(inputs.input_vcfs[0])"
                    ]
                }
            ], 
            "inputs": [
                {
                    "id": "#bcftools_isec.cwl/input_vcfs", 
                    "type": {
                        "type": "array", 
                        "items": "File"
                    }, 
                    "inputBinding": {
                        "position": 1
                    }
                }, 
                {
                    "id": "#bcftools_isec.cwl/collapse_records", 
                    "type": "string", 
                    "default": "all", 
                    "inputBinding": {
                        "prefix": "-c"
                    }
                }, 
                {
                    "id": "#bcftools_isec.cwl/nfiles", 
                    "type": "int", 
                    "default": null, 
                    "inputBinding": {
                        "prefix": "-n"
                    }
                }, 
                {
                    "id": "#bcftools_isec.cwl/output_files", 
                    "type": "int", 
                    "default": null, 
                    "inputBinding": {
                        "prefix": "-w"
                    }
                }, 
                {
                    "id": "#bcftools_isec.cwl/output_prefix", 
                    "type": "string", 
                    "default": "pref", 
                    "inputBinding": {
                        "prefix": "-p", 
                        "valueFrom": "$(self + '_' + inputs.input_vcfs[0].nameroot)"
                    }
                }, 
                {
                    "id": "#bcftools_isec.cwl/output_type", 
                    "type": "string", 
                    "default": "z", 
                    "inputBinding": {
                        "prefix": "-O"
                    }
                }
            ], 
            "outputs": [
                {
                    "id": "#bcftools_isec.cwl/isec_vcf", 
                    "type": "File", 
                    "outputBinding": {
                        "glob": "$(inputs.output_prefix + '_' + inputs.input_vcfs[0].nameroot + '/' + \"*.vcf.gz\")"
                    }
                }
            ], 
            "id": "#bcftools_isec.cwl"
        }, 
        {
            "class": "CommandLineTool", 
            "label": "remove indels from vcf", 
            "baseCommand": [
                "bcftools", 
                "isec"
            ], 
            "requirements": [
                {
                    "class": "DockerRequirement", 
                    "dockerPull": "quay.io/biocontainers/bcftools:1.8--1"
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "class": "InitialWorkDirRequirement", 
                    "listing": [
                        "$(inputs.input_vcf)"
                    ]
                }
            ], 
            "inputs": [
                {
                    "id": "#bcftools_isec_pairwise.cwl/input_vcf", 
                    "type": "File", 
                    "inputBinding": {
                        "position": 1
                    }
                }, 
                {
                    "id": "#bcftools_isec_pairwise.cwl/final_vcf", 
                    "type": "File", 
                    "inputBinding": {
                        "position": 2
                    }
                }, 
                {
                    "id": "#bcftools_isec_pairwise.cwl/collapse_records", 
                    "type": "string", 
                    "default": "all", 
                    "inputBinding": {
                        "prefix": "-c"
                    }
                }, 
                {
                    "id": "#bcftools_isec_pairwise.cwl/nfiles", 
                    "type": "int", 
                    "default": null, 
                    "inputBinding": {
                        "prefix": "-n"
                    }
                }, 
                {
                    "id": "#bcftools_isec_pairwise.cwl/output_files", 
                    "type": "int", 
                    "default": null, 
                    "inputBinding": {
                        "prefix": "-w"
                    }
                }, 
                {
                    "id": "#bcftools_isec_pairwise.cwl/output_prefix", 
                    "type": "string", 
                    "default": "pref", 
                    "inputBinding": {
                        "prefix": "-o", 
                        "valueFrom": "$(inputs.input_vcf.nameroot + '.vcf.gz')"
                    }
                }, 
                {
                    "id": "#bcftools_isec_pairwise.cwl/output_type", 
                    "type": "string", 
                    "default": "z", 
                    "inputBinding": {
                        "prefix": "-O"
                    }
                }
            ], 
            "outputs": [
                {
                    "id": "#bcftools_isec_pairwise.cwl/isec_vcf", 
                    "type": "File", 
                    "outputBinding": {
                        "glob": "$(inputs.input_vcf.nameroot + \"*.vcf.gz\")"
                    }
                }
            ], 
            "id": "#bcftools_isec_pairwise.cwl"
        }, 
        {
            "class": "CommandLineTool", 
            "label": "merge VCFs together", 
            "baseCommand": [
                "bcftools", 
                "merge"
            ], 
            "requirements": [
                {
                    "class": "DockerRequirement", 
                    "dockerPull": "quay.io/biocontainers/bcftools:1.8--1"
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "class": "InitialWorkDirRequirement", 
                    "listing": [
                        "$(inputs.input_vcfs[0])"
                    ]
                }
            ], 
            "inputs": [
                {
                    "id": "#bcftools_merge.cwl/input_vcfs", 
                    "type": {
                        "type": "array", 
                        "items": "File"
                    }, 
                    "inputBinding": {
                        "position": 1
                    }
                }, 
                {
                    "id": "#bcftools_merge.cwl/force-samples", 
                    "type": "boolean", 
                    "default": true, 
                    "inputBinding": {
                        "prefix": "--force-samples"
                    }
                }
            ], 
            "outputs": [
                {
                    "id": "#bcftools_merge.cwl/merged_vcf", 
                    "type": "stdout"
                }
            ], 
            "stdout": "merged.vcf", 
            "id": "#bcftools_merge.cwl"
        }, 
        {
            "class": "CommandLineTool", 
            "label": "remove indels from vcf", 
            "baseCommand": [
                "bcftools", 
                "reheader"
            ], 
            "requirements": [
                {
                    "class": "DockerRequirement", 
                    "dockerPull": "quay.io/biocontainers/bcftools:1.8--1"
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "class": "InitialWorkDirRequirement", 
                    "listing": [
                        "$(inputs.input_vcf)"
                    ]
                }
            ], 
            "inputs": [
                {
                    "id": "#bcftools_reheader.cwl/input_vcf", 
                    "type": "File", 
                    "inputBinding": {
                        "position": 1
                    }
                }, 
                {
                    "id": "#bcftools_reheader.cwl/samples", 
                    "type": "File", 
                    "inputBinding": {
                        "prefix": "-s"
                    }
                }, 
                {
                    "id": "#bcftools_reheader.cwl/output_name", 
                    "type": "string"
                }
            ], 
            "outputs": [
                {
                    "id": "#bcftools_reheader.cwl/reheaded_vcf", 
                    "type": "stdout"
                }
            ], 
            "stdout": "$(inputs.output_name)", 
            "id": "#bcftools_reheader.cwl"
        }, 
        {
            "class": "Workflow", 
            "label": "Merge and reheader VCF files", 
            "requirements": [
                {
                    "class": "StepInputExpressionRequirement"
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "class": "SubworkflowFeatureRequirement"
                }
            ], 
            "inputs": [
                {
                    "id": "#merge_reheader_workflow.cwl/input_vcfs", 
                    "type": {
                        "type": "array", 
                        "items": "File"
                    }
                }, 
                {
                    "id": "#merge_reheader_workflow.cwl/force_samples", 
                    "type": "boolean", 
                    "default": true
                }, 
                {
                    "id": "#merge_reheader_workflow.cwl/sample_file_name", 
                    "type": "string", 
                    "default": "samples.txt"
                }, 
                {
                    "id": "#merge_reheader_workflow.cwl/output_name", 
                    "type": "string", 
                    "default": "collapsed.vcf"
                }
            ], 
            "outputs": [
                {
                    "id": "#merge_reheader_workflow.cwl/merged_reheaded_vcf", 
                    "type": "File", 
                    "outputSource": "#merge_reheader_workflow.cwl/reheader_vcf/reheaded_vcf"
                }
            ], 
            "steps": [
                {
                    "id": "#merge_reheader_workflow.cwl/create_samples_file", 
                    "run": "#printf.cwl", 
                    "in": [
                        {
                            "id": "#merge_reheader_workflow.cwl/create_samples_file/contents", 
                            "source": "#merge_reheader_workflow.cwl/input_vcfs", 
                            "valueFrom": "${\n  var paths = []\n\n  for (var i = 0; i < self.length; i++) {\n      paths.push(self[i].basename);\n  }\n\n  console.log(paths);\n\n  return paths.join('\\n')\n}\n"
                        }, 
                        {
                            "id": "#merge_reheader_workflow.cwl/create_samples_file/file_name", 
                            "source": "#merge_reheader_workflow.cwl/sample_file_name"
                        }
                    ], 
                    "out": [
                        "#merge_reheader_workflow.cwl/create_samples_file/printed"
                    ]
                }, 
                {
                    "id": "#merge_reheader_workflow.cwl/merge_vcfs", 
                    "run": "#bcftools_merge.cwl", 
                    "in": [
                        {
                            "id": "#merge_reheader_workflow.cwl/merge_vcfs/input_vcfs", 
                            "source": "#merge_reheader_workflow.cwl/input_vcfs"
                        }, 
                        {
                            "id": "#merge_reheader_workflow.cwl/merge_vcfs/force_samples", 
                            "source": "#merge_reheader_workflow.cwl/force_samples"
                        }
                    ], 
                    "out": [
                        "#merge_reheader_workflow.cwl/merge_vcfs/merged_vcf"
                    ]
                }, 
                {
                    "id": "#merge_reheader_workflow.cwl/reheader_vcf", 
                    "run": "#bcftools_reheader.cwl", 
                    "in": [
                        {
                            "id": "#merge_reheader_workflow.cwl/reheader_vcf/input_vcf", 
                            "source": "#merge_reheader_workflow.cwl/merge_vcfs/merged_vcf"
                        }, 
                        {
                            "id": "#merge_reheader_workflow.cwl/reheader_vcf/samples", 
                            "source": "#merge_reheader_workflow.cwl/create_samples_file/printed"
                        }, 
                        {
                            "id": "#merge_reheader_workflow.cwl/reheader_vcf/output_name", 
                            "source": "#merge_reheader_workflow.cwl/output_name"
                        }
                    ], 
                    "out": [
                        "#merge_reheader_workflow.cwl/reheader_vcf/reheaded_vcf"
                    ]
                }
            ], 
            "id": "#merge_reheader_workflow.cwl"
        }, 
        {
            "class": "Workflow", 
            "label": "CompSNP", 
            "doc": "Comparative SNP calling pipeline which determines the best SNP candidates with\nwhich to discrimate between samples.\n", 
            "requirements": [
                {
                    "class": "SubworkflowFeatureRequirement"
                }, 
                {
                    "class": "StepInputExpressionRequirement"
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "class": "MultipleInputFeatureRequirement"
                }, 
                {
                    "class": "ScatterFeatureRequirement"
                }
            ], 
            "inputs": [
                {
                    "id": "#main/input_fastq", 
                    "type": {
                        "type": "array", 
                        "items": "File"
                    }
                }, 
                {
                    "id": "#main/reference", 
                    "type": "File"
                }, 
                {
                    "id": "#main/threads", 
                    "type": "int", 
                    "default": 1
                }, 
                {
                    "id": "#main/minimum_mapping_quality", 
                    "type": "int", 
                    "default": 20
                }, 
                {
                    "id": "#main/exclude_alignments_by_flag", 
                    "type": "int", 
                    "default": 2304
                }, 
                {
                    "id": "#main/disc_min_avg_bq", 
                    "type": "int", 
                    "default": 10
                }, 
                {
                    "id": "#main/disc_min_var_freq", 
                    "type": "float", 
                    "default": 0.5
                }, 
                {
                    "id": "#main/strict_min_avg_bq", 
                    "type": "int", 
                    "default": 10
                }, 
                {
                    "id": "#main/strict_min_var_freq", 
                    "type": "float", 
                    "default": 0.8
                }
            ], 
            "outputs": [
                {
                    "id": "#main/vcf", 
                    "type": "File", 
                    "outputSource": "#main/merge_reheader_strict_vcfs/merged_reheaded_vcf"
                }
            ], 
            "steps": [
                {
                    "id": "#main/minimap2_align", 
                    "doc": "First we align the input sample FASTQs to the supplied reference\n", 
                    "run": "#alignment_workflow.cwl", 
                    "in": [
                        {
                            "id": "#main/minimap2_align/input_reads", 
                            "source": "#main/input_fastq"
                        }, 
                        {
                            "id": "#main/minimap2_align/reference", 
                            "source": "#main/reference"
                        }, 
                        {
                            "id": "#main/minimap2_align/threads", 
                            "source": "#main/threads"
                        }
                    ], 
                    "out": [
                        "#main/minimap2_align/aligned_SAM"
                    ]
                }, 
                {
                    "id": "#main/sam_2_bam", 
                    "doc": "Filtering of the alignments performed in the previous step occur at this stage,\nalong with conversion to BAM\n", 
                    "run": "#sam2bam_workflow.cwl", 
                    "scatter": "#main/sam_2_bam/input_sam", 
                    "in": [
                        {
                            "id": "#main/sam_2_bam/input_sam", 
                            "source": "#main/minimap2_align/aligned_SAM"
                        }, 
                        {
                            "id": "#main/sam_2_bam/threads", 
                            "source": "#main/threads"
                        }, 
                        {
                            "id": "#main/sam_2_bam/mapping_quality", 
                            "source": "#main/minimum_mapping_quality"
                        }, 
                        {
                            "id": "#main/sam_2_bam/exclude_alignments_by_flag", 
                            "source": "#main/exclude_alignments_by_flag"
                        }
                    ], 
                    "out": [
                        "#main/sam_2_bam/sorted_BAM"
                    ]
                }, 
                {
                    "id": "#main/varscan2_discovery", 
                    "doc": "In this step, we called SNPs with a relatively low threshold for allele frequency\nin order to collect the set of common SNP candidates.\n", 
                    "run": "#varscan2.cwl", 
                    "scatter": "#main/varscan2_discovery/input_bam", 
                    "in": [
                        {
                            "id": "#main/varscan2_discovery/input_bam", 
                            "source": "#main/sam_2_bam/sorted_BAM"
                        }, 
                        {
                            "id": "#main/varscan2_discovery/out_name", 
                            "valueFrom": "$('discovery_varscan2.vcf')"
                        }, 
                        {
                            "id": "#main/varscan2_discovery/input_ref", 
                            "source": "#main/reference"
                        }, 
                        {
                            "id": "#main/varscan2_discovery/min_avg_bq", 
                            "source": "#main/disc_min_avg_bq"
                        }, 
                        {
                            "id": "#main/varscan2_discovery/min_var_freq", 
                            "source": "#main/disc_min_var_freq"
                        }
                    ], 
                    "out": [
                        "#main/varscan2_discovery/compressed_indexed_vcf"
                    ]
                }, 
                {
                    "id": "#main/isec_discovery", 
                    "doc": "In this step, we collapse our discovery SNP calling VCF files so that we have\na list of all common SNP candidates in one file.\n", 
                    "run": "#bcftools_isec.cwl", 
                    "in": [
                        {
                            "id": "#main/isec_discovery/input_vcfs", 
                            "source": "#main/varscan2_discovery/compressed_indexed_vcf"
                        }, 
                        {
                            "id": "#main/isec_discovery/nfiles", 
                            "source": "#main/varscan2_discovery/compressed_indexed_vcf", 
                            "valueFrom": "$(self.length)"
                        }, 
                        {
                            "id": "#main/isec_discovery/output_files", 
                            "valueFrom": "$(1)"
                        }
                    ], 
                    "out": [
                        "#main/isec_discovery/isec_vcf"
                    ]
                }, 
                {
                    "id": "#main/index_discovery_vcf", 
                    "doc": "Lets quickly re-index the recently made discovery_collapsed_vcf and move on.\n", 
                    "run": "#bcftools_index.cwl", 
                    "in": [
                        {
                            "id": "#main/index_discovery_vcf/input_vcf", 
                            "source": "#main/isec_discovery/isec_vcf"
                        }
                    ], 
                    "out": [
                        "#main/index_discovery_vcf/indexed_vcf"
                    ]
                }, 
                {
                    "id": "#main/varscan2_strict", 
                    "doc": "Now we re-call the SNPs with a higher stringency in terms of allele frequency,\nin order to get a list, for each sample, of the likely SNPs relative to the reference.\n", 
                    "run": "#varscan2.cwl", 
                    "scatter": "#main/varscan2_strict/input_bam", 
                    "in": [
                        {
                            "id": "#main/varscan2_strict/input_bam", 
                            "source": "#main/sam_2_bam/sorted_BAM"
                        }, 
                        {
                            "id": "#main/varscan2_strict/out_name", 
                            "valueFrom": "$('strict_varscan2.vcf')"
                        }, 
                        {
                            "id": "#main/varscan2_strict/input_ref", 
                            "source": "#main/reference"
                        }, 
                        {
                            "id": "#main/varscan2_strict/min_avg_bq", 
                            "source": "#main/strict_min_avg_bq"
                        }, 
                        {
                            "id": "#main/varscan2_strict/min_var_freq", 
                            "source": "#main/strict_min_var_freq"
                        }
                    ], 
                    "out": [
                        "#main/varscan2_strict/compressed_indexed_vcf"
                    ]
                }, 
                {
                    "id": "#main/isec_strict_discovery", 
                    "doc": "Next, for each sample we find any SNPs found during strict SNP calling which do not\nappear in the discovery SNP candidates list.\n", 
                    "run": "#bcftools_isec_pairwise.cwl", 
                    "scatter": "#main/isec_strict_discovery/input_vcf", 
                    "in": [
                        {
                            "id": "#main/isec_strict_discovery/input_vcf", 
                            "source": "#main/varscan2_strict/compressed_indexed_vcf"
                        }, 
                        {
                            "id": "#main/isec_strict_discovery/final_vcf", 
                            "source": "#main/index_discovery_vcf/indexed_vcf"
                        }, 
                        {
                            "id": "#main/isec_strict_discovery/nfiles", 
                            "valueFrom": "$(1)"
                        }, 
                        {
                            "id": "#main/isec_strict_discovery/output_files", 
                            "valueFrom": "$(1)"
                        }
                    ], 
                    "out": [
                        "#main/isec_strict_discovery/isec_vcf"
                    ]
                }, 
                {
                    "id": "#main/index_strict_vcfs", 
                    "doc": "Lets quickly re-index the recently made strict isec_vcfs and move on.\n", 
                    "run": "#bcftools_index.cwl", 
                    "scatter": "#main/index_strict_vcfs/input_vcf", 
                    "in": [
                        {
                            "id": "#main/index_strict_vcfs/input_vcf", 
                            "source": "#main/isec_strict_discovery/isec_vcf"
                        }
                    ], 
                    "out": [
                        "#main/index_strict_vcfs/indexed_vcf"
                    ]
                }, 
                {
                    "id": "#main/merge_reheader_strict_vcfs", 
                    "doc": "Finally we merge the strict VCFs and reheader them so that each sample can be distinguished.\n", 
                    "run": "#merge_reheader_workflow.cwl", 
                    "in": [
                        {
                            "id": "#main/merge_reheader_strict_vcfs/input_vcfs", 
                            "source": "#main/index_strict_vcfs/indexed_vcf"
                        }
                    ], 
                    "out": [
                        "#main/merge_reheader_strict_vcfs/merged_reheaded_vcf"
                    ]
                }
            ], 
            "id": "#main"
        }, 
        {
            "class": "CommandLineTool", 
            "label": "compress a file using bgzip from htslib", 
            "baseCommand": [
                "bgzip"
            ], 
            "requirements": [
                {
                    "class": "DockerRequirement", 
                    "dockerPull": "quay.io/biocontainers/htslib:1.7--0"
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "class": "InitialWorkDirRequirement", 
                    "listing": [
                        "$(inputs.input)"
                    ]
                }
            ], 
            "inputs": [
                {
                    "id": "#htslib_bgzip.cwl/input", 
                    "type": "File", 
                    "inputBinding": {
                        "position": 1
                    }
                }
            ], 
            "outputs": [
                {
                    "id": "#htslib_bgzip.cwl/compressed", 
                    "type": "File", 
                    "outputBinding": {
                        "glob": "*.vcf.gz"
                    }
                }
            ], 
            "id": "#htslib_bgzip.cwl"
        }, 
        {
            "class": "Workflow", 
            "label": "Alignment Evaluation Pipeline - minimap2", 
            "doc": "Align the input against the reference genome with minimap2", 
            "requirements": [
                {
                    "class": "ScatterFeatureRequirement"
                }
            ], 
            "inputs": [
                {
                    "id": "#alignment_workflow.cwl/input_reads", 
                    "type": {
                        "type": "array", 
                        "items": "File"
                    }
                }, 
                {
                    "id": "#alignment_workflow.cwl/reference", 
                    "type": "File"
                }, 
                {
                    "id": "#alignment_workflow.cwl/threads", 
                    "type": "int", 
                    "default": 1
                }
            ], 
            "outputs": [
                {
                    "id": "#alignment_workflow.cwl/aligned_SAM", 
                    "type": {
                        "type": "array", 
                        "items": "File"
                    }, 
                    "outputSource": "#alignment_workflow.cwl/align/SAM"
                }
            ], 
            "steps": [
                {
                    "id": "#alignment_workflow.cwl/minimap2_index", 
                    "run": "#minimap2_index.cwl", 
                    "in": [
                        {
                            "id": "#alignment_workflow.cwl/minimap2_index/reference", 
                            "source": "#alignment_workflow.cwl/reference"
                        }
                    ], 
                    "out": [
                        "#alignment_workflow.cwl/minimap2_index/indexes"
                    ]
                }, 
                {
                    "id": "#alignment_workflow.cwl/align", 
                    "run": "#minimap2.cwl", 
                    "scatter": "#alignment_workflow.cwl/align/input_reads", 
                    "in": [
                        {
                            "id": "#alignment_workflow.cwl/align/input_reads", 
                            "source": "#alignment_workflow.cwl/input_reads"
                        }, 
                        {
                            "id": "#alignment_workflow.cwl/align/reference", 
                            "source": "#alignment_workflow.cwl/minimap2_index/indexes"
                        }, 
                        {
                            "id": "#alignment_workflow.cwl/align/threads", 
                            "source": "#alignment_workflow.cwl/threads"
                        }
                    ], 
                    "out": [
                        "#alignment_workflow.cwl/align/SAM"
                    ]
                }
            ], 
            "id": "#alignment_workflow.cwl"
        }, 
        {
            "class": "CommandLineTool", 
            "label": "minimap2", 
            "doc": "Align input against reference with minimap2", 
            "baseCommand": "minimap2", 
            "requirements": [
                {
                    "class": "DockerRequirement", 
                    "dockerPull": "quay.io/biocontainers/minimap2:2.10--1"
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "class": "ShellCommandRequirement"
                }
            ], 
            "inputs": [
                {
                    "id": "#minimap2.cwl/output_sam", 
                    "type": "boolean", 
                    "default": true, 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "-a"
                    }
                }, 
                {
                    "id": "#minimap2.cwl/input_reads", 
                    "type": "File", 
                    "inputBinding": {
                        "position": 3
                    }
                }, 
                {
                    "id": "#minimap2.cwl/reference", 
                    "type": "File", 
                    "inputBinding": {
                        "position": 2
                    }
                }, 
                {
                    "id": "#minimap2.cwl/preset", 
                    "type": "string", 
                    "default": "map-ont", 
                    "inputBinding": {
                        "prefix": "-x"
                    }
                }, 
                {
                    "id": "#minimap2.cwl/threads", 
                    "type": "int", 
                    "default": 1, 
                    "inputBinding": {
                        "prefix": "-t"
                    }
                }
            ], 
            "outputs": [
                {
                    "id": "#minimap2.cwl/SAM", 
                    "type": "File", 
                    "outputBinding": {
                        "glob": "*.sam"
                    }
                }
            ], 
            "stdout": "$(inputs.input_reads.nameroot + '_minimap2' + '.sam')", 
            "id": "#minimap2.cwl"
        }, 
        {
            "class": "CommandLineTool", 
            "label": "minimap2 index", 
            "doc": "Create index files reference genome with minimap2", 
            "baseCommand": [
                "minimap2"
            ], 
            "requirements": [
                {
                    "class": "DockerRequirement", 
                    "dockerPull": "quay.io/biocontainers/minimap2:2.10--1"
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "class": "InitialWorkDirRequirement", 
                    "listing": [
                        "$(inputs.reference)"
                    ]
                }
            ], 
            "inputs": [
                {
                    "id": "#minimap2_index.cwl/reference", 
                    "type": "File", 
                    "inputBinding": {
                        "position": 1, 
                        "valueFrom": "$(self.basename)"
                    }
                }
            ], 
            "arguments": [
                {
                    "prefix": "-d", 
                    "valueFrom": "$(inputs.reference.basename + \".minimap2.idx\")"
                }
            ], 
            "outputs": [
                {
                    "id": "#minimap2_index.cwl/indexes", 
                    "type": "File", 
                    "outputBinding": {
                        "glob": "$(inputs.reference.basename + \".minimap2.idx\")"
                    }
                }
            ], 
            "id": "#minimap2_index.cwl"
        }, 
        {
            "class": "Workflow", 
            "label": "SAM2BAM", 
            "doc": "Convert SAM to a sorted and indexed BAM", 
            "requirements": [
                {
                    "class": "MultipleInputFeatureRequirement"
                }, 
                {
                    "class": "StepInputExpressionRequirement"
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "class": "SubworkflowFeatureRequirement"
                }
            ], 
            "inputs": [
                {
                    "id": "#sam2bam_workflow.cwl/input_sam", 
                    "type": "File"
                }, 
                {
                    "id": "#sam2bam_workflow.cwl/threads", 
                    "type": "int"
                }, 
                {
                    "id": "#sam2bam_workflow.cwl/mapping_quality", 
                    "type": "int", 
                    "default": 0
                }, 
                {
                    "id": "#sam2bam_workflow.cwl/exclude_alignments_by_flag", 
                    "type": "int", 
                    "default": null
                }
            ], 
            "outputs": [
                {
                    "id": "#sam2bam_workflow.cwl/sorted_BAM", 
                    "type": "File", 
                    "outputSource": "#sam2bam_workflow.cwl/sort/sorted_BAM"
                }
            ], 
            "steps": [
                {
                    "id": "#sam2bam_workflow.cwl/sam2bam", 
                    "run": "#samtools_view.cwl", 
                    "in": [
                        {
                            "id": "#sam2bam_workflow.cwl/sam2bam/input_sam", 
                            "source": "#sam2bam_workflow.cwl/input_sam"
                        }, 
                        {
                            "id": "#sam2bam_workflow.cwl/sam2bam/threads", 
                            "source": "#sam2bam_workflow.cwl/threads"
                        }, 
                        {
                            "id": "#sam2bam_workflow.cwl/sam2bam/mapping_quality", 
                            "source": "#sam2bam_workflow.cwl/mapping_quality"
                        }, 
                        {
                            "id": "#sam2bam_workflow.cwl/sam2bam/exclude_alignments_by_flag", 
                            "source": "#sam2bam_workflow.cwl/exclude_alignments_by_flag"
                        }
                    ], 
                    "out": [
                        "#sam2bam_workflow.cwl/sam2bam/output_bam"
                    ]
                }, 
                {
                    "id": "#sam2bam_workflow.cwl/sort", 
                    "run": "#samtools_sort.cwl", 
                    "in": [
                        {
                            "id": "#sam2bam_workflow.cwl/sort/input_bam", 
                            "source": "#sam2bam_workflow.cwl/sam2bam/output_bam"
                        }, 
                        {
                            "id": "#sam2bam_workflow.cwl/sort/output_bam_path", 
                            "source": "#sam2bam_workflow.cwl/input_sam", 
                            "valueFrom": "$(self.nameroot + '.bam')"
                        }, 
                        {
                            "id": "#sam2bam_workflow.cwl/sort/threads", 
                            "source": "#sam2bam_workflow.cwl/threads"
                        }
                    ], 
                    "out": [
                        "#sam2bam_workflow.cwl/sort/sorted_BAM"
                    ]
                }
            ], 
            "id": "#sam2bam_workflow.cwl"
        }, 
        {
            "class": "CommandLineTool", 
            "label": "Call SNPs using VarScan2", 
            "baseCommand": [
                "samtools", 
                "mpileup"
            ], 
            "requirements": [
                {
                    "class": "DockerRequirement", 
                    "dockerPull": "quay.io/biocontainers/samtools:latest"
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "class": "InitialWorkDirRequirement", 
                    "listing": [
                        "$(inputs.input_fa)"
                    ]
                }
            ], 
            "inputs": [
                {
                    "id": "#samtools_mpileup.cwl/input_fa", 
                    "type": "File", 
                    "inputBinding": {
                        "prefix": "-f"
                    }
                }, 
                {
                    "id": "#samtools_mpileup.cwl/input_bam", 
                    "type": "File", 
                    "inputBinding": {
                        "position": 1
                    }
                }, 
                {
                    "id": "#samtools_mpileup.cwl/min-bq", 
                    "type": "int", 
                    "default": 10, 
                    "inputBinding": {
                        "prefix": "-Q"
                    }
                }, 
                {
                    "id": "#samtools_mpileup.cwl/min-mq", 
                    "type": "int", 
                    "default": 20, 
                    "inputBinding": {
                        "prefix": "-q"
                    }
                }, 
                {
                    "id": "#samtools_mpileup.cwl/no_baq", 
                    "type": "boolean", 
                    "default": true, 
                    "inputBinding": {
                        "prefix": "-B"
                    }
                }, 
                {
                    "id": "#samtools_mpileup.cwl/region", 
                    "type": "string", 
                    "default": null, 
                    "inputBinding": {
                        "prefix": "--region"
                    }
                }
            ], 
            "outputs": [
                {
                    "id": "#samtools_mpileup.cwl/pileup", 
                    "streamable": true, 
                    "type": "stdout"
                }
            ], 
            "stdout": "$(inputs.input_bam.nameroot)", 
            "id": "#samtools_mpileup.cwl"
        }, 
        {
            "class": "CommandLineTool", 
            "label": "SAMtools sort", 
            "doc": "Sorts a BAM file", 
            "baseCommand": [
                "samtools", 
                "sort"
            ], 
            "requirements": [
                {
                    "class": "DockerRequirement", 
                    "dockerPull": "quay.io/biocontainers/samtools:latest"
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }
            ], 
            "inputs": [
                {
                    "id": "#samtools_sort.cwl/input_bam", 
                    "type": "File", 
                    "inputBinding": {
                        "position": 1
                    }
                }, 
                {
                    "id": "#samtools_sort.cwl/output_bam_path", 
                    "type": "string", 
                    "inputBinding": {
                        "prefix": "-o"
                    }
                }, 
                {
                    "id": "#samtools_sort.cwl/threads", 
                    "type": "int", 
                    "default": 1, 
                    "inputBinding": {
                        "prefix": "-@"
                    }
                }, 
                {
                    "id": "#samtools_sort.cwl/by_readname", 
                    "type": "boolean", 
                    "default": false, 
                    "inputBinding": {
                        "prefix": "-n"
                    }
                }
            ], 
            "outputs": [
                {
                    "id": "#samtools_sort.cwl/sorted_BAM", 
                    "type": "File", 
                    "outputBinding": {
                        "glob": "*.bam"
                    }
                }
            ], 
            "id": "#samtools_sort.cwl"
        }, 
        {
            "class": "CommandLineTool", 
            "label": "SAMtools view", 
            "doc": "Converts a SAM to BAM", 
            "baseCommand": [
                "samtools", 
                "view"
            ], 
            "requirements": [
                {
                    "class": "DockerRequirement", 
                    "dockerPull": "quay.io/biocontainers/samtools:latest"
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }
            ], 
            "inputs": [
                {
                    "id": "#samtools_view.cwl/input_sam", 
                    "type": "File", 
                    "inputBinding": {
                        "position": 1
                    }
                }, 
                {
                    "id": "#samtools_view.cwl/output_bam_path", 
                    "type": "string", 
                    "default": "tmp.bam", 
                    "inputBinding": {
                        "prefix": "-o"
                    }
                }, 
                {
                    "id": "#samtools_view.cwl/threads", 
                    "type": "int", 
                    "default": 1, 
                    "inputBinding": {
                        "prefix": "-@"
                    }
                }, 
                {
                    "id": "#samtools_view.cwl/include_header_flag", 
                    "type": "boolean", 
                    "default": true, 
                    "inputBinding": {
                        "prefix": "-h"
                    }
                }, 
                {
                    "id": "#samtools_view.cwl/output_bam_flag", 
                    "type": "boolean", 
                    "default": true, 
                    "inputBinding": {
                        "prefix": "-b"
                    }
                }, 
                {
                    "id": "#samtools_view.cwl/mapping_quality", 
                    "type": "int", 
                    "default": 0, 
                    "inputBinding": {
                        "prefix": "-q"
                    }
                }, 
                {
                    "id": "#samtools_view.cwl/exclude_alignments_by_flag", 
                    "type": "int", 
                    "default": null, 
                    "inputBinding": {
                        "prefix": "-F"
                    }
                }
            ], 
            "outputs": [
                {
                    "id": "#samtools_view.cwl/output_bam", 
                    "type": "File", 
                    "outputBinding": {
                        "glob": "*.bam"
                    }
                }
            ], 
            "id": "#samtools_view.cwl"
        }, 
        {
            "class": "CommandLineTool", 
            "doc": "printf to file", 
            "baseCommand": [
                "printf"
            ], 
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ], 
            "inputs": [
                {
                    "id": "#printf.cwl/contents", 
                    "type": "string", 
                    "inputBinding": {
                        "position": 1
                    }
                }, 
                {
                    "id": "#printf.cwl/file_name", 
                    "type": "string"
                }
            ], 
            "outputs": [
                {
                    "id": "#printf.cwl/printed", 
                    "type": "stdout"
                }
            ], 
            "stdout": "$(inputs.file_name)", 
            "id": "#printf.cwl"
        }, 
        {
            "class": "CommandLineTool", 
            "label": "Call SNPs using VarScan2", 
            "baseCommand": [
                "varscan", 
                "mpileup2snp"
            ], 
            "requirements": [
                {
                    "class": "DockerRequirement", 
                    "dockerPull": "quay.io/biocontainers/varscan:2.4.3--0"
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }
            ], 
            "inputs": [
                {
                    "id": "#mpileup2snp.cwl/mpileup", 
                    "type": "File", 
                    "streamable": true
                }, 
                {
                    "id": "#mpileup2snp.cwl/name", 
                    "type": "string", 
                    "default": "varscan2.vcf"
                }, 
                {
                    "id": "#mpileup2snp.cwl/min_avg_bq", 
                    "type": "int", 
                    "default": 15, 
                    "inputBinding": {
                        "prefix": "--min-avg-qual"
                    }
                }, 
                {
                    "id": "#mpileup2snp.cwl/strand_filter", 
                    "type": "boolean", 
                    "default": true, 
                    "inputBinding": {
                        "prefix": "--strand-filter", 
                        "valueFrom": "$(self | 0)"
                    }
                }, 
                {
                    "id": "#mpileup2snp.cwl/min_var_freq", 
                    "type": "float", 
                    "default": 0.01, 
                    "inputBinding": {
                        "prefix": "--min-var-freq"
                    }
                }, 
                {
                    "id": "#mpileup2snp.cwl/variants_only", 
                    "type": "boolean", 
                    "default": false, 
                    "inputBinding": {
                        "prefix": "--variants", 
                        "valueFrom": "$(self | 0)"
                    }
                }, 
                {
                    "id": "#mpileup2snp.cwl/min_p_val", 
                    "type": "float", 
                    "default": null, 
                    "inputBinding": {
                        "prefix": "--p-value"
                    }
                }
            ], 
            "arguments": [
                "--output-vcf 1"
            ], 
            "stdin": "$(inputs.mpileup.path)", 
            "outputs": [
                {
                    "id": "#mpileup2snp.cwl/vcf", 
                    "type": "stdout"
                }
            ], 
            "stdout": "$(inputs.name)", 
            "id": "#mpileup2snp.cwl"
        }, 
        {
            "class": "Workflow", 
            "label": "Call SNPs using VarScan2", 
            "requirements": [
                {
                    "class": "StepInputExpressionRequirement"
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "class": "SubworkflowFeatureRequirement"
                }
            ], 
            "inputs": [
                {
                    "id": "#varscan2.cwl/input_bam", 
                    "type": "File"
                }, 
                {
                    "id": "#varscan2.cwl/input_ref", 
                    "type": "File"
                }, 
                {
                    "id": "#varscan2.cwl/out_name", 
                    "type": "string", 
                    "default": "varscan2.vcf"
                }, 
                {
                    "id": "#varscan2.cwl/min_avg_bq", 
                    "type": "int", 
                    "default": 6
                }, 
                {
                    "id": "#varscan2.cwl/strand_filter", 
                    "type": "boolean", 
                    "default": true
                }, 
                {
                    "id": "#varscan2.cwl/min_var_freq", 
                    "type": "float", 
                    "default": 0.25
                }, 
                {
                    "id": "#varscan2.cwl/variants_only", 
                    "type": "boolean", 
                    "default": false
                }, 
                {
                    "id": "#varscan2.cwl/min_p_val", 
                    "type": "float", 
                    "default": null
                }
            ], 
            "outputs": [
                {
                    "id": "#varscan2.cwl/compressed_indexed_vcf", 
                    "type": "File", 
                    "outputSource": "#varscan2.cwl/index_vcf/indexed_vcf"
                }
            ], 
            "steps": [
                {
                    "id": "#varscan2.cwl/mpileup", 
                    "run": "#samtools_mpileup.cwl", 
                    "in": [
                        {
                            "id": "#varscan2.cwl/mpileup/input_bam", 
                            "source": "#varscan2.cwl/input_bam"
                        }, 
                        {
                            "id": "#varscan2.cwl/mpileup/input_fa", 
                            "source": "#varscan2.cwl/input_ref"
                        }
                    ], 
                    "out": [
                        "#varscan2.cwl/mpileup/pileup"
                    ]
                }, 
                {
                    "id": "#varscan2.cwl/mpileup2snp", 
                    "run": "#mpileup2snp.cwl", 
                    "in": [
                        {
                            "id": "#varscan2.cwl/mpileup2snp/mpileup", 
                            "source": "#varscan2.cwl/mpileup/pileup"
                        }, 
                        {
                            "id": "#varscan2.cwl/mpileup2snp/name", 
                            "source": "#varscan2.cwl/out_name", 
                            "valueFrom": "$(inputs.mpileup.nameroot + '_' + self)"
                        }, 
                        {
                            "id": "#varscan2.cwl/mpileup2snp/min_avg_bq", 
                            "source": "#varscan2.cwl/min_avg_bq"
                        }, 
                        {
                            "id": "#varscan2.cwl/mpileup2snp/strand_filter", 
                            "source": "#varscan2.cwl/strand_filter"
                        }, 
                        {
                            "id": "#varscan2.cwl/mpileup2snp/min_var_freq", 
                            "source": "#varscan2.cwl/min_var_freq"
                        }, 
                        {
                            "id": "#varscan2.cwl/mpileup2snp/variants_only", 
                            "source": "#varscan2.cwl/variants_only"
                        }, 
                        {
                            "id": "#varscan2.cwl/mpileup2snp/min_p_val", 
                            "source": "#varscan2.cwl/min_p_val"
                        }
                    ], 
                    "out": [
                        "#varscan2.cwl/mpileup2snp/vcf"
                    ]
                }, 
                {
                    "id": "#varscan2.cwl/compress_vcf", 
                    "run": "#htslib_bgzip.cwl", 
                    "in": [
                        {
                            "id": "#varscan2.cwl/compress_vcf/input", 
                            "source": "#varscan2.cwl/mpileup2snp/vcf"
                        }
                    ], 
                    "out": [
                        "#varscan2.cwl/compress_vcf/compressed"
                    ]
                }, 
                {
                    "id": "#varscan2.cwl/index_vcf", 
                    "run": "#bcftools_index.cwl", 
                    "in": [
                        {
                            "id": "#varscan2.cwl/index_vcf/input_vcf", 
                            "source": "#varscan2.cwl/compress_vcf/compressed"
                        }
                    ], 
                    "out": [
                        "#varscan2.cwl/index_vcf/indexed_vcf"
                    ]
                }
            ], 
            "id": "#varscan2.cwl"
        }
    ]
}