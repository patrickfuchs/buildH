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


help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
.PHONY: help

