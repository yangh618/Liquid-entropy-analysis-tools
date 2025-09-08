#!/bin/bash

# For this run, POSCAR, POTCAR, XDATCAR are required.

#  Pre-run
# ============================================
cat XDATCAR | pdfxdat_cubic 1> /tmp/Nul
pdf2s2_v2 -x XDATCAR > tot.s1

