#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/epitopeprediction
========================================================================================
    Github : https://github.com/nf-core/epitopeprediction
    Website: https://nf-co.re/epitopeprediction
    Slack  : https://nfcore.slack.com/channels/epitopeprediction
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { EPITOPEPREDICTION } from './workflows/epitopeprediction'

//
// WORKFLOW: Run main nf-core/epitopeprediction analysis pipeline
//
workflow NFCORE_EPITOPEPREDICTION {
    EPITOPEPREDICTION ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_EPITOPEPREDICTION ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
