default: help

install-env: ## Install conda environement
	conda env create -f binder/environment.yml
.PHONY: install-env


update-env: init ## Update conda environement
	conda env update -f binder/environment-dev.yml
.PHONY: update-env


tests: ## Run tests
	pytest tests
.PHONY: tests


lint: ## Lint code
	pycodestyle buildh \
	&& pydocstyle buildh \
	&& pylint buildh
.PHONY: lint


docs: ## Build documentation
	cd docs && make html
.PHONY: docs


# see doc and tutorial:
# https://packaging.python.org/tutorials/packaging-projects/
build: ## Build package for PyPI
	python -m build
.PHONY: build


upload-to-pypi: ## Upload to PyPI
	python -m twine upload dist/*
	# clean packages
	rm -f dist/*.tar.gz dist/*.whl dist/*.egg
.PHONY: upload-to-pypi


help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
.PHONY: help

