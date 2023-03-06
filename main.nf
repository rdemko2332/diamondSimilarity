#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//---------------------------------------------------------------
// Param Checking 
//---------------------------------------------------------------

if(params.seqFile) {
  seqs = Channel.fromPath( params.seqFile ).splitFasta( by:params.fastaSubsetSize, file:true  )
}
else {
  throw new Exception("Missing params.seqFile")
}

if (params.preConfiguredDatabase) {
  if (!params.database) {
    throw new Exception("Missing params.database")
  }
}

//--------------------------------------------------------------------------
// Includes
//--------------------------------------------------------------------------

include { nonConfiguredDatabase } from './modules/diamondSimilarity.nf'
include { preConfiguredDatabase } from './modules/diamondSimilarity.nf'

//--------------------------------------------------------------------------
// Main Workflow
//--------------------------------------------------------------------------

workflow {
  
  if (!params.preConfiguredDatabase) {
    nonConfiguredDatabase(seqs)
  }
   
  else if (params.preConfiguredDatabase) {
    preConfiguredDatabase(seqs)
  }

}

