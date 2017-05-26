test:
	pytest -v --cov=nuclmm --cov-report term-missing nuclmm/*.py nuclmm/tests/*.py

style:
	pep8 nuclmm/*.py nuclmm/tests/*.py

loc:
	cloc --exclude-list-file=<(echo nuclmm/_version.py) nuclmm/*.py
	cloc nuclmm/tests/test_*.py

clean:
	rm -rf __pycache__/
