# Variables
SHELL = /bin/bash -O globstar
.SHELLFLAGS = -o pipefail -c
CPUS = 10
MEMS = 100
MEM_PER_CPU = $(shell echo '$(MEMS)/$(CPUS)' | bc)

# Programs
SLURM = srun --chdir=$(abspath .) --cpus-per-task=$(CPUS) --export=ALL --job-name="$(PROJECT)_$(@F)" --mem=$(MEMS)G --mail-user="jaewoong@localhost" --mail-type=END,FAIL
PYTHON = /home/Live/jaewoong/Lineage_tracing-Onboarding/.venv/bin/python3 -B -u -X faulthandler
R = /usr/bin/env RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE /usr/bin/Rscript
rwildcard = $(foreach d,$(wildcard $(1:=/*)),$(call rwildcard,$d,$2) $(filter $(subst *,%,$2),$d))
product = $(foreach x,$1,$(addprefix $x.,$2))

# References
DARLIN_DIR = /home/Store/MethylTree/metadata/DARLIN_HSC/downstream_R
DARLIN_SAMPLE = $(DARLIN_DIR)/df_sample/sample_sheet.tsv.gz
DARLIN_H5AD = $(addprefix $(DARLIN_DIR)/out_dir/mouse_,$(addsuffix _neutral_thresh_DARLIN_HSC_BM_0.66_0.92_readcutoff3_minCov0.1_met_ExcludeSex0.h5ad,$(MICE)))

# Constants
PROJECT = Lineage
MICE = LL731 LL901 LL653E1 LL653E6
OPTIMAL = Darlin
DATASETS = Test Darlin
SEX_MODES = IncludeSex ExcludeSex

# Options
EXCLUDE_IncludeSex =
EXCLUDE_ExcludeSex = chrX chrY chrM
LOGGING = 1>$(abspath $@).log 2>&1
LOGGING_DIR = 1>$(abspath $(@D)).log 2>&1
LOGGING_TEE = 2>&1 | tee $(abspath $@).log

all: test optimal step01 step02 step03# step04
	@date | mail -s "[$(PROJECT)] Making $@ has been finished on $(shell hostname)" "jaewoong@gist.ac.kr"
.PHONY: all

clean:
.PHONY: clean

test: step01-test step02-test step03-test# step04-test
	@date | mail -s "[$(PROJECT)] Making $@ has been finished on $(shell hostname)" "jaewoong@gist.ac.kr"
.PHONY: test

optimal: step01-optimal step02-optimal step03-optimal# step04-optimal
	@date | mail -s "[$(PROJECT)] Making $@ has been finished on $(shell hostname)" "jaewoong@gist.ac.kr"
.PHONY: optimal

# Step 01 (Scanpy import)
Output/01_Scanpy_import:
	mkdir -p $@

Output/01_Scanpy_import/Darlin.h5ad: Python/import_h5ad.py $(DARLIN_H5AD) $(DARLIN_SAMPLE) | Output/01_Scanpy_import
	$(SLURM) $(PYTHON) $(abspath $^ $@) --names $(MICE) --clone "large_clone_id" --join "outer" $(LOGGING)

Output/01_Scanpy_import/Test.h5ad: Python/import_h5ad.py /home/Store/MethylTree/metadata/DARLIN_HSC/downstream_R/out_dir/mouse_LL653E1_neutral_thresh_DARLIN_HSC_BM_0.66_0.92_readcutoff3_minCov0.1_met_ExcludeSex0.h5ad $(DARLIN_SAMPLE) | Output/01_Scanpy_import
	$(SLURM) $(PYTHON) $(abspath $^ $@) --names LL653E1 --clone "large_clone_id" --join "outer" $(LOGGING)

step01-test: Output/01_Scanpy_import/Test.h5ad
step01-optimal: Output/01_Scanpy_import/$(OPTIMAL).h5ad
step01: $(addprefix Output/01_Scanpy_import/,$(addsuffix .h5ad,$(DATASETS)))
.PHONY: step01-test step01-optimal step01

# Step 02 (QC covered regions)
Output/02_QC_covered_regions:
	mkdir -p $@

Output/02_QC_covered_regions/Figures:
	mkdir -p $@

Output/02_QC_covered_regions/%.cell.tsv.gz: Python/qc_covered_regions.py Output/01_Scanpy_import/%.h5ad | Output/02_QC_covered_regions
	$(SLURM) $(PYTHON) $(abspath $^ $@) --target "cell" $(LOGGING)

Output/02_QC_covered_regions/%.region.tsv.gz: Python/qc_covered_regions.py Output/01_Scanpy_import/%.h5ad | Output/02_QC_covered_regions
	$(SLURM) $(PYTHON) $(abspath $^ $@) --target "region" $(LOGGING)

Output/02_QC_covered_regions/Figures/%.zip: Python/draw_covered_regions.py Output/02_QC_covered_regions/%.cell.tsv.gz Output/02_QC_covered_regions/%.region.tsv.gz | Output/02_QC_covered_regions/Figures
	$(SLURM) $(PYTHON) $(abspath $^ $@) $(LOGGING)

step02-test: Output/02_QC_covered_regions/Test.cell.tsv.gz Output/02_QC_covered_regions/Test.region.tsv.gz Output/02_QC_covered_regions/Figures/Test.zip
step02-optimal: Output/02_QC_covered_regions/$(OPTIMAL).cell.tsv.gz Output/02_QC_covered_regions/$(OPTIMAL).region.tsv.gz Output/02_QC_covered_regions/Figures/$(OPTIMAL).zip
step02: $(addprefix Output/02_QC_covered_regions/,$(addsuffix .cell.tsv.gz,$(DATASETS))) $(addprefix Output/02_QC_covered_regions/,$(addsuffix .region.tsv.gz,$(DATASETS))) $(addprefix Output/02_QC_covered_regions/Figures/,$(addsuffix .zip,$(DATASETS)))
.PHONY: step02-test step02-optimal step02

# Step 03 (MethylTree per individual)
Output/03_MethylTree:
	mkdir -p $@

Output/03_MethylTree/Figures:
	mkdir -p $@

Output/03_MethylTree/%.IncludeSex.h5ad: Python/run_methyltree.py Output/01_Scanpy_import/%.h5ad | Output/03_MethylTree
	$(SLURM) $(PYTHON) $(abspath $^ $@) --clone "large_clone_id" --exclude $(EXCLUDE_IncludeSex) $(LOGGING)

Output/03_MethylTree/%.ExcludeSex.h5ad: Python/run_methyltree.py Output/01_Scanpy_import/%.h5ad | Output/03_MethylTree
	$(SLURM) $(PYTHON) $(abspath $^ $@) --clone "large_clone_id" --exclude $(EXCLUDE_ExcludeSex) $(LOGGING)

Output/03_MethylTree/%.tsv.gz: Python/summarize_methyltree.py $(addprefix Output/03_MethylTree/%.,$(addsuffix .h5ad,$(SEX_MODES))) | Output/03_MethylTree
	$(SLURM) $(PYTHON) $(abspath $^ $@) $(LOGGING)

Output/03_MethylTree/Figures/%.zip: Python/draw_methyltree.py $(addprefix Output/03_MethylTree/%.,$(addsuffix .h5ad,$(SEX_MODES))) | Output/03_MethylTree/Figures
	$(SLURM) $(PYTHON) $(abspath $^ $@) $(LOGGING)

step03-test: $(addprefix Output/03_MethylTree/Test.,$(addsuffix .h5ad,$(SEX_MODES))) Output/03_MethylTree/Test.tsv.gz Output/03_MethylTree/Figures/Test.zip
step03-optimal: $(addprefix Output/03_MethylTree/$(OPTIMAL).,$(addsuffix .h5ad,$(SEX_MODES))) Output/03_MethylTree/$(OPTIMAL).tsv.gz Output/03_MethylTree/Figures/$(OPTIMAL).zip
step03: $(call product,$(addprefix Output/03_MethylTree/,$(DATASETS)),$(addsuffix .h5ad,$(SEX_MODES))) $(addprefix Output/03_MethylTree/,$(addsuffix .tsv.gz,$(DATASETS))) $(addprefix Output/03_MethylTree/Figures/,$(addsuffix .zip,$(DATASETS)))
.PHONY: step03-test step03-optimal step03

# Step 04 (Coverage down-sampling)
Output/04_Downsampling:
	mkdir -p $@

Output/04_Downsampling/Figures:
	mkdir -p $@

Output/04_Downsampling/Darlin.%.tsv.gz: Python/downsample_methyltree.py Output/01_Scanpy_import/Darlin.h5ad | Output/04_Downsampling
	$(SLURM) $(PYTHON) $(abspath $^ $@) $(DOWNSAMPLING_OPTIONS) $(LOGGING)

Output/04_Downsampling/Test.%.tsv.gz: Python/downsample_methyltree.py Output/01_Scanpy_import/Test.h5ad | Output/04_Downsampling
	$(SLURM) $(PYTHON) $(abspath $^ $@) $(DOWNSAMPLING_OPTIONS) $(LOGGING)

Output/04_Downsampling/Figures/Darlin.zip: Python/draw_downsampling.py $(addprefix Output/04_Downsampling/Darlin.,$(addsuffix .tsv.gz,$(call product,$(MICE),$(KEEP_PERCENTS)))) | Output/04_Downsampling/Figures
	$(SLURM) $(PYTHON) $(abspath $^ $@) $(LOGGING)

Output/04_Downsampling/Figures/Test.zip: Python/draw_downsampling.py $(addprefix Output/04_Downsampling/Test.LL653E1.,$(addsuffix .tsv.gz,$(TEST_KEEP_PERCENTS))) | Output/04_Downsampling/Figures
	$(SLURM) $(PYTHON) $(abspath $^ $@) $(LOGGING)

step04-test: Output/04_Downsampling/Figures/Test.zip
step04-optimal: Output/04_Downsampling/Figures/$(OPTIMAL).zip
step04: $(addprefix Output/04_Downsampling/Figures/,$(addsuffix .zip,$(DATASETS)))
.PHONY: step04-test step04-optimal step04
