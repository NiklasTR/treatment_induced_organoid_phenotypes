.PHONY: clean data lint requirements sync_data_to_s3 sync_data_from_s3

#################################################################################
# GLOBALS                                                                       #
#################################################################################

PROJECT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
BUCKET = [OPTIONAL] your-bucket-for-syncing-data (do not include 's3://')
PROFILE = default
PROJECT_NAME = promise
PYTHON_INTERPRETER = python3

ifeq (,$(shell which conda))
HAS_CONDA=False
else
HAS_CONDA=True
endif

#################################################################################
# COMMANDS                                                                      #
#################################################################################

## Pull Dependencies
# TODO fill out dependicy pulling
dependencies:
	cd src/data
	git clone https://github.com/NiklasTR/phenotypespectrum.git\
	cd ../models
	git clone https://github.com/NiklasTR/phenotypespectrum.git
# also docker dependency managed via singularity singularity pull docker://gtca/mofa2	

## Install Python Dependencies
requirements: test_environment
	$(PYTHON_INTERPRETER) -m pip install -U pip setuptools wheel
	$(PYTHON_INTERPRETER) -m pip install -r requirements.txt

## Make Dataset
data: requirements
	$(PYTHON_INTERPRETER) src/data/make_dataset.py data/raw data/processed

## Delete all compiled Python files
clean:
	find . -type f -name "*.py[co]" -delete
	find . -type d -name "__pycache__" -delete

## Lint using flake8
lint:
	flake8 src

## Upload Data to S3
sync_data_to_s3:
ifeq (default,$(PROFILE))
	aws s3 sync data/ s3://$(BUCKET)/data/
else
	aws s3 sync data/ s3://$(BUCKET)/data/ --profile $(PROFILE)
endif

## Download Data from S3
sync_data_from_s3:
ifeq (default,$(PROFILE))
	aws s3 sync s3://$(BUCKET)/data/ data/
else
	aws s3 sync s3://$(BUCKET)/data/ data/ --profile $(PROFILE)
endif

## Set up python interpreter environment
create_environment:
ifeq (True,$(HAS_CONDA))
		@echo ">>> Detected conda, creating conda environment."
ifeq (3,$(findstring 3,$(PYTHON_INTERPRETER)))
	conda create --name $(PROJECT_NAME) python=3
else
	conda create --name $(PROJECT_NAME) python=2.7
endif
		@echo ">>> New conda env created. Activate with:\nsource activate $(PROJECT_NAME)"
else
	$(PYTHON_INTERPRETER) -m pip install -q virtualenv virtualenvwrapper
	@echo ">>> Installing virtualenvwrapper if not already installed.\nMake sure the following lines are in shell startup file\n\
	export WORKON_HOME=$$HOME/.virtualenvs\nexport PROJECT_HOME=$$HOME/Devel\nsource /usr/local/bin/virtualenvwrapper.sh\n"
	@bash -c "source `which virtualenvwrapper.sh`;mkvirtualenv $(PROJECT_NAME) --python=$(PYTHON_INTERPRETER)"
	@echo ">>> New virtualenv created. Activate with:\nworkon $(PROJECT_NAME)"
endif

## Test python environment is setup correctly
test_environment: dependencies
	$(PYTHON_INTERPRETER) test_environment.py

## Set up python interpreter environment
export_environment:
ifeq (True,$(HAS_CONDA))
		@echo ">>> Detected conda, exporting conda environment."
		conda list --explicit > conda_env.txt
endif
		@echo ">>> conda env exported. Import with:\nconda env create --file conda_env.txt --name $(PROJECT_NAME)"


#################################################################################
# CLUSTER RULES                                                                 #
#################################################################################

## This Makefile is located on the group's cluster (data4 = "cold storage") and has can be executed with make sync 
#sync_all:
#        # current hot storage location of the promise project
#        rsync -rav /data3/promise /data4/b110_image_archive/HC1092/promise_coldstorage
#        # legacy sources of code and data
#        rsync -rav /data/promise /data4/b110_image_archive/HC1092/archive/b110_cluster_promise
#        rsync -rav /data2/valentint/promise /data4/b110_image_archive/HC1092/archive/b110_cluster_valentini
#        rsync -rav /data3/promise_legacy /data4/b110_image_archive/HC1092/archive/isilon2
#        # additional legacy sources of code and data that are not easily accessed programatically
#        # NOT RUN rsync -rav rindtorf@odcf-lsf01.dkfz.de: rindtorf@b110-sc2cn01:/data4/b110_image_archive/HC1092/archive/odcf_rstudio
#        # rsync -rav Documents/GitHub/promise rindtorf@b110-sc2cn01:/data4/b110_image_archive/HC1092/archive/rindtorff_macbook

#sync:
#        rsync -rav /data3/promise /data4/b110_image_archive/HC1092/promise_coldstorage

#sync_legacy:
#        rsync -rav /data/promise /data4/b110_image_archive/HC1092/archive/b110_cluster_promise
#        rsync -rav /data2/valentint/promise /data4/b110_image_archive/HC1092/archive/b110_cluster_valentini  
#        rsync -rav /data3/promise_legacy /data4/b110_image_archive/HC1092/archive/isilon2

# TODO: remove if not needed any longer
# rsync -rav rindtorf@b110-sc2cn01:/data/valentini/promise/PROMISE/data-10x-4t-c-16z/hdf5projection/ data/raw/PROMISE/data-10x-4t-c-16z/hdf5projection/

# TODO needs implementation of one-line ssh and sync
# backup:
#	ssh rindtorf@b110-sc2cn01; cd /data4/b110_image_archive/HC1092; make sync

### access
ssh_b110:
	ssh rindtorf@b110-sc2cn01

ssh_bsub: .log
	ssh rindtorf@odcf-lsf01.dkfz.de
	cd /dkfz/groups/shared/OE0049/B110-Isilon2/promise/
	touch .log
	
### modules
load_modules:
	module load libpng/1.6.37
	module load gdal/3.0.2
	module load pandoc/2.2.1; 
	#module load anaconda3/2019.07
	module load R/4.0.0 
	module load gcc/7.2.0 
	module load anaconda3/2020.11
	module load tmux/3.1c

# ssh user@socket command < /path/to/file/on/local/machine

#################################################################################
# PROJECT RULES                                                                 #
#################################################################################



## Run Live-Dead-Classifier

## Run Line Difference Analysis

## Run DrugEffects Analysis

### filter dead organoids from all organoids in processedFeatures

#################################################################################
# Self Documenting Commands                                                     #
#################################################################################

.DEFAULT_GOAL := help

# Inspired by <http://marmelab.com/blog/2016/02/29/auto-documented-makefile.html>
# sed script explained:
# /^##/:
# 	* save line in hold space
# 	* purge line
# 	* Loop:
# 		* append newline + line to hold space
# 		* go to next line
# 		* if line starts with doc comment, strip comment character off and loop
# 	* remove target prerequisites
# 	* append hold space (+ newline) to line
# 	* replace newline plus comments by `---`
# 	* print line
# Separate expressions are necessary because labels cannot be delimited by
# semicolon; see <http://stackoverflow.com/a/11799865/1968>
.PHONY: help
help:
	@echo "$$(tput bold)Available rules:$$(tput sgr0)"
	@echo
	@sed -n -e "/^## / { \
		h; \
		s/.*//; \
		:doc" \
		-e "H; \
		n; \
		s/^## //; \
		t doc" \
		-e "s/:.*//; \
		G; \
		s/\\n## /---/; \
		s/\\n/ /g; \
		p; \
	}" ${MAKEFILE_LIST} \
	| LC_ALL='C' sort --ignore-case \
	| awk -F '---' \
		-v ncol=$$(tput cols) \
		-v indent=19 \
		-v col_on="$$(tput setaf 6)" \
		-v col_off="$$(tput sgr0)" \
	'{ \
		printf "%s%*s%s ", col_on, -indent, $$1, col_off; \
		n = split($$2, words, " "); \
		line_length = ncol - indent; \
		for (i = 1; i <= n; i++) { \
			line_length -= length(words[i]) + 1; \
			if (line_length <= 0) { \
				line_length = ncol - indent - length(words[i]) - 1; \
				printf "\n%*s ", -indent, " "; \
			} \
			printf "%s ", words[i]; \
		} \
		printf "\n"; \
	}' \
	| more $(shell test $(shell uname) = Darwin && echo '--no-init --raw-control-chars')
