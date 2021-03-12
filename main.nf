// Hybrid assembly workflow: DSL2

nextflow.enable.dsl=2

// Helper functions

def get_single_fastx( glob ){
    return channel.fromPath(glob) | map { file -> tuple(file.baseName, file) }
}

def get_paired_fastq( glob ){
    return channel.fromFilePairs(glob, flat: true)
}

def get_matching_data( channel1, channel2, illumina = false ){

    // Get matching data by ID (first field) from two channels
    // by crossing, checking for matched ID and returning
    // a combined data tuple with a single ID (first field)

    channel1.cross(channel2).map { crossed ->
        if (crossed[0][0] == crossed[1][0]){
            if (illumina){
                // First channel returning Illumina only for hybrid reference assembly
                crossed[0]
            } else {
                // Return mix of channels for evaluations in hybrid assembly workflow
                tuple( crossed[0][0], *crossed[0][1..-1], *crossed[1][1..-1] )
            }
        } else {
            null
        }
    }
    .filter { it != null }

}

// ONT Quality control subworkflow


params.workflow = 'hybrid'
params.reference = ""
params.fastq = 'fastq/*.fq'
params.outdir = 'results'
params.length = 200
params.quality = 7
params.coverage = 200
params.genome_size = '2.8m'
params.preset = 'minimap2-ont'  // coverm

// ONT assembly and polishing subworkflow

params.assembly_options = "--plasmids"
params.medaka_model = "r941_min_high_g360"

if ( file(params.medaka_model).exists() ){
    params.medaka_model = file(params.medaka_model)
}

// Illumina assembly, hybrid correction and reference comparison subworkflow

params.fasta = "fasta/*.fasta"
params.illumina = "fastq/*_R{1,2}.fastq.gz"
params.depth = 200
params.assembler = "skesa"

params.tag = null
params.saureus = true
params.kpneumoniae = false

// Stage files

reference = file(params.reference)

// Modules 

include { Nanoq } from './modules/nanoq'
include { Rasusa } from './modules/rasusa'
include { CoverM  } from './modules/coverm'
include { NanoqStatistics } from './modules/nanoq'
include { Fastp } from './modules/fastp'
include { Shovill } from './modules/shovill'
include { Flye } from './modules/flye'
include { Racon } from './modules/racon'
include { Medaka } from './modules/medaka'
inlcude { UnicyclerHybrid } from './modules/unicycler'

// Assign tags for output names / folders and parsing with NanoPath

include { Pilon as PilonCorrection } from './modules/pilon' addParams( tag: 'medaka' )
include { Dnadiff as ONTComparison } from './modules/dnadiff' addParams( tag: 'medaka' )
include { Dnadiff as HybridComparison } from './modules/dnadiff' addParams( tag: 'medaka_hybrid' )
include { Genotype as IlluminaGenotype } from './modules/genotype' addParams( tag: 'illumina' )
include { Genotype as AssemblyGenotype } from './modules/genotype' addParams( tag: 'preassembled' )
include { Genotype as HybridGenotype } from './modules/genotype' addParams( tag: 'hybrid' )
include { Genotype as MedakaGenotype } from './modules/genotype' addParams( tag: 'ont' )



workflow ont_qc {
    take:
        reads // id, reads
    main:
        NanoqStatistics(reads)
        Nanoq(reads)
        CoverM(Nanoq.out, reference)
        Rasusa (Nanoq.out) 
    emit:
        Rasusa.out // id, reads
}


workflow ont_assembly {
    take:
        reads // id, reads
    main:
        Flye(reads)
        get_matching_data(Flye.out, reads, false) | Racon
        Medaka(Racon.out)
        MedakaGenotype(Medaka.out)
    emit:
        Flye.out  // id, assembly
        Medaka.out  // id, polished assembly
}


workflow illumina_assembly {
    take:
        reads // id, fwd, rev
    main:
        Fastp(reads) | Shovill | IlluminaGenotype
    emit:
        Fastp.out  // id, fwd, rev
        Shovill.out  // id, fasta
}

workflow hybrid_correction {
    take:
        illumina_reads // id, fwd, rev
        ont_assembly // id, fasta
        reference_assembly // id, fasta
    main:
        get_matching_data(ont_assembly, illumina_reads, false) | PilonCorrection
        get_matching_data(ont_assembly, reference_assembly, false) | ONTComparison
        get_matching_data(PilonCorrection.out, reference_assembly, false) | HybridComparison
        HybridGenotype(PilonCorrection.out)  // Genotype corrected assembly
    emit:
        PilonCorrection.out
}       


workflow np_core_assembly {

   if (params.workflow == "genotype") {
       // Genotyping from pre-assembled genomes only
       get_single_fastx(params.fasta)  | AssemblyGenotype
   }  
   if (params.workflow == "ont"){
        // ONT standard workflow with Flye + Medaka and genotyping
        get_single_fastx(params.fastq) | ont_qc | ont_assembly
   }
   if (params.workflow == "illumina") {
       // Illumina standard workflow and genotyping
       get_paired_fastq(params.illumina) | illumina_assembly
   }
   if (params.workflow == "hybrid-pilon") {
        // ONT and Illumina reference assemblies
        get_single_fastx(params.fastq) | ont_qc | ont_assembly
        get_matching_data(get_paired_fastq(params.illumina), get_single_fastx(params.fastq), true) | illumina_assembly
        
        hybrid_correction(
            illumina_assembly.out[0], // qc reads
            ont_assembly.out[1], // polished ont assembly
            illumina_assembly.out[1] // reference illumina assembly
        )
   }
   if (params.workflow == "hybrid-unicycler") {
        // Unicycler long-read hybrid assembly
        get_single_fastx(params.fastq) | ont_qc
        get_paired_fastq(params.illumina) | Fastp
        get_matching_data(Fastp.out, ont_qc.out, false) | UnicyclerHybrid
        
   }
}

// Execute

workflow {
    np_core_assembly()
}