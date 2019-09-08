configfile: "config.yaml"

import glob
import os
import logging

SAMPLES, = glob_wildcards("data/{sample}.rds")
print("=====")
print(SAMPLES)
print("=====")

rule all:
	input:
	  expand("output/{sample}.html", sample = SAMPLES),
	  expand("output/settings_{sample}.rds", sample = SAMPLES),
	  expand("output/vp_{sample}.rds", sample = SAMPLES),
	  expand("output/rel_diff_abund_{sample}.rds", sample = SAMPLES),
	  expand("output/cms_{sample}.rds", sample = SAMPLES),
	  expand("output/cms_summary_{sample}.rds", sample = SAMPLES),
	  expand("output/res_{sample}.rds", sample = SAMPLES),
	  expand("output/de_summary_{sample}.rds", sample = SAMPLES),
	  expand("output/sim_{sample}.rds", sample = SAMPLES),
	  expand("output/sim_summary_{sample}.rds", sample = SAMPLES)

## -------------------------------------------------------------------------- ##
## Reports
## -------------------------------------------------------------------------- ##


rule Rmarkdown:
	input:
		"data/{sample}.rds",
		"output/settings_{sample}.rds",
		"analysis/report.Rmd"
	output:
		"output/{sample}.html"
	script:
		"analysis/report.Rmd"


rule settings:
    input:
        "data/{sample}.rds"
    output:
        "output/settings_{sample}.rds"
    shell:
        "Rscript --vanilla analysis/settings.R {input} {output}"

rule variancePartition:
    input:
        sce="data/{sample}.rds", settings = "output/settings_{sample}.rds"
    output:
        vp="output/vp_{sample}.rds"
    shell:
        "Rscript --vanilla analysis/variancePartition.R {input.sce} {input.settings} {output.vp}"


rule rel_abund_diff:
    input:
        sce="data/{sample}.rds", settings="output/settings_{sample}.rds"
    output:
        rel_diff_abund="output/rel_diff_abund_{sample}.rds"
    shell:
        "Rscript --vanilla analysis/rel_diff_abund.R {input.sce} {input.settings} {output.rel_diff_abund}"

rule cms:
    input:
        sce="data/{sample}.rds", settings="output/settings_{sample}.rds"
    output:
        cms="output/cms_{sample}.rds", cms_summary="output/cms_summary_{sample}.rds"
    shell:
        "Rscript --vanilla analysis/cms.R {input.sce} {input.settings} {output.cms} {output.cms_summary}"

rule de_analysis:
    input:
        sce="data/{sample}.rds", settings="output/settings_{sample}.rds"
    output:
        res="output/res_{sample}.rds", de_summary="output/de_summary_{sample}.rds"
    shell:
        "Rscript --vanilla analysis/de_analysis.R {input.sce} {input.settings} {output.res} {output.de_summary}"

rule simulation:
    input:
        res="output/res_{sample}.rds", de_summary="output/de_summary_{sample}.rds", sce="data/{sample}.rds", summary_vp="output/vp_{sample}.rds", settings="output/settings_{sample}.rds"
    output:
        sim="output/sim_{sample}.rds"
    shell:
        "Rscript --vanilla analysis/simulation.R {input.res} {input.de_summary} {input.sce} {input.summary_vp} {input.settings} {output.sim}"


rule summary:
    input:
        summary_vp="output/vp_{sample}.rds", rel_diff_abund="output/rel_diff_abund_{sample}.rds", de_summary="output/de_summary_{sample}.rds", cms_summary="output/cms_summary_{sample}.rds", sim="output/sim_{sample}.rds", sce="data/{sample}.rds", settings="output/settings_{sample}.rds"
    output:
        sim_summary="output/sim_summary_{sample}.rds"
    shell:
        "Rscript --vanilla analysis/summary.R {input.summary_vp} {input.rel_diff_abund} {input.de_summary} {input.cms_summary} {input.sim} {input.sce} {input.settings} {output.sim_summary}"
