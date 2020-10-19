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
params.fastq = 'fastq/*.fq'
params.outdir = 'results'
params.reference = "$PWD/jdk.fasta"
params.length = 200
params.quality = 7
params.coverage = 200
params.genome_size = '2.8m'
params.preset = 'minimap2-ont'  // coverm



include { Nanoq } from './modules/nanoq'
include { Rasusa } from './modules/rasusa'
include { CoverM  } from './modules/coverm'
include { NanoqStatistics } from './modules/nanoq'

workflow ont_qc {
    take:
        reads // id, reads
    main:
        NanoqStatistics(reads)
        Nanoq(reads)
        CoverM(Nanoq.out, params.reference)
        Rasusa (Nanoq.out) 
    emit:
        Rasusa.out // id, reads
}

// ONT assembly and polishing subworkflow

params.assembly_options = "--plasmids"
params.medaka_model = "r941_min_high_g360"

include { Flye } from './modules/flye'
include { Racon } from './modules/racon'
include { Medaka } from './modules/medaka'

workflow ont_assembly {
    take:
        reads // id, reads
    main:
        Flye(reads)
        get_matching_data(Flye.out, reads) | Racon
        Medaka(Racon.out)
        MedakaGenotype(Medaka.out)
    emit:
        Flye.out  // id, assembly
        Medaka.out  // id, polished assembly
}

// Illumina assembly, hybrid correction and reference comparison subworkflow

params.fasta = "fasta/*.fasta"
params.illumina = "fastq/*_R{1,2}.fastq.gz"
params.depth = 200
params.assembler = "skesa"

params.tag = null
params.saureus = true
params.kpneumoniae = false

include { Fastp } from './modules/fastp'
include { Shovill } from './modules/shovill'

// Assign tags for separate output folders in params.outdir / hybrid / params.tag
include { Pilon as FlyePilon } from './modules/pilon' addParams( tag: 'flye' )
include { Pilon as MedakaPilon } from './modules/pilon' addParams( tag: 'medaka' )

include { Dnadiff as MedakaComparison } from './modules/dnadiff' addParams( tag: 'medaka' )
include { Dnadiff as FlyeHybridComparison } from './modules/dnadiff' addParams( tag: 'flye_hybrid' )
include { Dnadiff as MedakaHybridComparison } from './modules/dnadiff' addParams( tag: 'medaka_hybrid' )

include { Genotype as IlluminaGenotype } from './modules/genotype' addParams( tag: 'illumina' )
include { Genotype as AssemblyGenotype } from './modules/genotype' addParams( tag: 'fasta' )
include { Genotype as MedakaHybridGenotype } from './modules/genotype' addParams( tag: 'hybrid' )
include { Genotype as MedakaGenotype } from './modules/genotype' addParams( tag: 'ont' )

workflow illumina_assembly {
    take:
        reads // id, fwd, rev
    main:
        Fastp(reads) | Shovill | IlluminaGenotype
    emit:
        Fastp.out
        Shovill.out
        IlluminaGenotype.out
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
   if (params.workflow == "hybrid") {
        // ONT and Illumina reference assemblies
        get_single_fastx(params.fastq) | ont_qc | ont_assembly
        
        get_matching_data(
            get_paired_fastq(params.illumina), get_single_fastx(params.fastq), true  // matching illumina only
        ) | illumina_assembly
        
        // Branch into hybrid corrections with Pilon
        get_matching_data(ont_assembly.out[0], illumina_assembly.out[0]) | FlyePilon  // flye assembly correction
        get_matching_data(ont_assembly.out[1], illumina_assembly.out[0]) | MedakaPilon  // polished assembly correction
        
        // TODO: Implement HyPo, POLCA and NextPolish alternatives, and chained polishing

        // Genotype corrected hybrid assemblies from Medaka 
        MedakaHybridGenotype(MedakaPilon.out)

        // Compare Flye and Medaka assemblies to Illumina reference assembly (DNADIFF)
        get_matching_data(ont_assembly.out[1], illumina_assembly.out[1]) | MedakaComparison // illumina vs. polished assembly
        get_matching_data(FlyePilon.out, illumina_assembly.out[1]) | FlyeHybridComparison // illumina vs. corrected raw assembly
        get_matching_data(MedakaPilon.out, illumina_assembly.out[1]) | MedakaHybridComparison // illumina vs. corrected polished assembly
   }
}

// Execute

workflow {
    np_core_assembly()
}