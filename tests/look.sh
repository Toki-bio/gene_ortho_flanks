#!/bin/bash
cd /home/toki/gof/toy
for g in FAHD2A_5140 PAM16_52 SQLE_137 UNC50_92; do echo "-- $g"; grep "^$g|[LR]|1000" out_keep/flanks.paf | cut -f1-12; grep "^$g|" out_keep/flanks.bed; done
