# Variables
SHELL = /bin/bash -O globstar
.SHELLFLAGS = -o pipefail -c
CPUS = 10
MEMS = 100
MEM_PER_CPU = $(shell echo '$(MEMS)/$(CPUS)' | bc)

# Programs
SLURM = srun --chdir=$(abspath .) --cpus-per-task=$(CPUS) --export=ALL --job-name="Lineage_$(@F)" --mem=$(MEMS)G --mail-user="jaewoong@gist.ac.kr" --mail-type=FAIL
PYTHON = /home/Live/jaewoong/Lineage_tracing-Onboarding/.venv/bin/python3 -B
R = /usr/bin/env RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE /usr/bin/Rscript
rwildcard = $(foreach d,$(wildcard $(1:=/*)),$(call rwildcard,$d,$2) $(filter $(subst *,%,$2),$d))
product = $(foreach x,$1,$(addprefix $x.,$2))

# Options
LOGGING = 1>$(abspath $@).log 2>&1
LOGGING_DIR = 1>$(abspath $(@D)).log 2>&1

all:
.PHONY: all

clean:
.PHONY: clean

test:
.PHONY: test

optimal:
.PHONY: optimal

# Step 01 ()
