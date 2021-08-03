tests:
	@echo "Running Unittests using tox"
	tox

dist:
	@echo "Creating Distributions"
        python setup.py sdist

deploy:
	@echo "Deploying package"
        twine upload dist/*

black:
	black --line-length 89 --exclude '(ursgal/wrapper_template/|.tox)' .
