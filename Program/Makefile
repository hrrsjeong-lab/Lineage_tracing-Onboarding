# Variables
SHELL = /bin/bash -O globstar
.SHELLFLAGS = -o pipefail -c
CPUS = 10
MEMS = 100
MEM_PER_CPU = $(shell echo '$(MEMS)/$(CPUS)' | bc)

# Programs
SLURM = srun --chdir=$(abspath .) --cpus-per-task=$(CPUS) --export=ALL --job-name="Lineage_$(@F)" --mem=$(MEMS)G --mail-user="jaewoong@gist.ac.kr" --mail-type=END,FAIL
PYTHON = /home/Live/jaewoong/Lineage_tracing-Onboarding/.venv/bin/python3 -B
R = /usr/bin/env RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE /usr/bin/Rscript
rwildcard = $(foreach d,$(wildcard $(1:=/*)),$(call rwildcard,$d,$2) $(filter $(subst *,%,$2),$d))
product = $(foreach x,$1,$(addprefix $x.,$2))

# References
DARLIN_DIR = /home/Store/MethylTree/metadata/DARLIN_HSC/downstream_R
DARLIN_SAMPLE = $(DARLIN_DIR)/df_sample/sample_sheet.tsv.gz
DARLIN_H5AD = $(addprefix $(DARLIN_DIR)/out_dir/mouse_,$(addsuffix _neutral_thresh_DARLIN_HSC_BM_0.66_0.92_readcutoff3_minCov0.1_met_ExcludeSex0.h5ad,$(MICE)))

# Constants
MICE = LL731 LL901 LL653E1 LL653E6
OPTIMAL = Darlin
DATASETS = Test Darlin

# Options
LOGGING = 1>$(abspath $@).log 2>&1
LOGGING_DIR = 1>$(abspath $(@D)).log 2>&1

all: test optimal step01
	@date | mail -s "[Slurm] Making $@ has been finished on $(shell hostname)" "jaewoong@gist.ac.kr"
.PHONY: all

clean:
.PHONY: clean

test: step01-test
	@date | mail -s "[Slurm] Making $@ has been finished on $(shell hostname)" "jaewoong@gist.ac.kr"
.PHONY: test

optimal: step01-optimal
	@date | mail -s "[Slurm] Making $@ has been finished on $(shell hostname)" "jaewoong@gist.ac.kr"
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
step01: $(addprefix Output/01_Scanpy_import/,$(addsuffix .h5ad,$(MICE)))
.PHONY: step01-test step01-optimal step01
