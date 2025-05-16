#!/usr/bin/bash
# multiqc on the mapping results

multiqc -n 043_r_multiqc_mouseMT_mapped_raw.html -f --title mapped_raw 042_d_STAR_map_raw/
