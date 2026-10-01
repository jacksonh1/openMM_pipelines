#################################################################################
# GLOBALS                                                                       #
#################################################################################

PROJECT_NAME = openmm_pipelines
ENV_NAME = openmm_env
# PYTHON_VERSION = 3.12
PYTHON_INTERPRETER = python

#################################################################################
# COMMANDS                                                                      #
#################################################################################


## Install/update conda environment
.PHONY: requirements
requirements:
	mamba env update --name $(ENV_NAME) --file environment.yml


## Delete all compiled Python files ("*.py[co]", "__pycache__", etc.)
.PHONY: clean
clean:
	find . -type f -name "*.py[co]" -delete
	find . -type d -name "__pycache__" -delete


## create conda environment
.PHONY: create_environment
create_environment:
	mamba env create --name $(ENV_NAME) -f environment.yml

	@echo ">>> conda env created. Activate with:\nconda activate $(ENV_NAME)"


#################################################################################
# PROJECT RULES                                                                 #
#################################################################################


## Format source code with black
.PHONY: format
format:
	black --config pyproject.toml $(PROJECT_NAME)


## Run tests
.PHONY: test
test:
	pytest tests


#################################################################################
# Self Documenting Commands                                                     #
#################################################################################

.DEFAULT_GOAL := help

define PRINT_HELP_PYSCRIPT
import re, sys; \
lines = '\n'.join([line for line in sys.stdin]); \
matches = re.findall(r'\n## (.*)\n[\s\S]+?\n([a-zA-Z_-]+):', lines); \
print('Available rules:\n'); \
print('\n'.join(['{:25}{}'.format(*reversed(match)) for match in matches]))
endef
export PRINT_HELP_PYSCRIPT

help:
	@$(PYTHON_INTERPRETER) -c "${PRINT_HELP_PYSCRIPT}" < $(MAKEFILE_LIST)
