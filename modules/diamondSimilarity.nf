#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process createDatabase {
  input:
    path newdbfasta

  output:
    path 'newdb.dmnd'

  script:
    template 'createDatabase.bash'
}

process diamondSimilarity {
  input:
    path fasta
    path database
    val pValCutoff 
    val lengthCutoff 
    val percentCutoff 
    val blastProgram  
    val blastArgs 
    
  output:
    path 'diamondSimilarity.out', emit: output_file
    path 'diamondSimilarity.log', emit: log_file

  script:
    template 'diamondSimilarity.bash'
}


workflow nonConfiguredDatabase {
  take:
    seqs

  main:
    database = createDatabase(params.databaseFasta)
    diamondSimilarityResults = diamondSimilarity(seqs, database, params.pValCutoff, params.lengthCutoff, params.percentCutoff, params.blastProgram, params.blastArgs)
    diamondSimilarityResults.output_file | collectFile(storeDir: params.outputDir, name: params.dataFile)
    diamondSimilarityResults.log_file | collectFile(storeDir: params.outputDir, name: params.logFile)
}

workflow preConfiguredDatabase {
  take:
    seqs

  main:
    diamondSimilarityResults = diamondSimilarity(seqs, params.database, params.pValCutoff, params.lengthCutoff, params.percentCutoff, params.blastProgram, params.blastArgs)
    diamondSimilarityResults.output_file | collectFile(storeDir: params.outputDir, name: params.dataFile)
    diamondSimilarityResults.log_file | collectFile(storeDir: params.outputDir, name: params.logFile)
    
}