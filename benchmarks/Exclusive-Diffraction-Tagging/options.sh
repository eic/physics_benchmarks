#!/bin/bash

## =============================================================================
## Shared step-control functions for Exclusive-Diffraction-Tagging benchmarks.
##
## Source this file (do not run it directly). It defines:
##
##   print_step_options_help  - prints the step-option help text; call this from
##                              the sourcing script's print_the_help function
##   parse_step_options "$@"  - parses step-control flags and sets:
##                                DO_ALL, DATA_INIT, DO_SIM, DO_REC, DO_ANALYSIS
##                              Calls print_the_help (must be defined by caller)
##                              on -h/--help or unknown options.
##
## Usage in a benchmark script:
##
##   source benchmarks/Exclusive-Diffraction-Tagging/options.sh
##
##   function print_the_help {
##     echo "USAGE: ${0} [--sim] [--rec] [--analysis] [--all]"
##     print_step_options_help
##     exit
##   }
##
##   parse_step_options "$@"
##
## =============================================================================

function print_step_options_help {
  echo "    The default behavior is to run all steps (sim, rec, analysis)."
  echo "STEP OPTIONS:"
  echo "  --data-init     Download the input event data"
  echo "  --sim, -s       Run the Geant4 simulation"
  echo "  --rec, -r       Run the reconstruction"
  echo "  --analysis, -a  Run the analysis scripts"
  echo "  --all           (default) Run all steps"
  echo "  -h, --help      Print this message"
}

function parse_step_options {
  DO_ALL=1
  DATA_INIT=
  DO_SIM=
  DO_REC=
  DO_ANALYSIS=

  while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
      -h|--help)
        print_the_help
        ;;
      --all)
        DO_ALL=1
        if [[ -n "${DO_REC}${DO_SIM}${DO_ANALYSIS}" ]]; then
          echo "Error: cannot use --all with other step arguments." 1>&2
          print_the_help
          exit 1
        fi
        shift
        ;;
      -s|--sim)
        DO_SIM=1
        DO_ALL=
        shift
        ;;
      --data-init)
        DATA_INIT=1
        DO_ALL=
        shift
        ;;
      -r|--rec)
        DO_REC=1
        DO_ALL=
        shift
        ;;
      -a|--analysis)
        DO_ANALYSIS=1
        DO_ALL=
        shift
        ;;
      *)
        echo "unknown option '$1'" 1>&2
        print_the_help
        exit 1
        ;;
    esac
  done
}
