test:
	pytest -v --cov=nuclmm --cov-report term-missing nuclmm/*.py nuclmm/tests/*.py

style:
	pep8 nuclmm/*.py nuclmm/tests/*.py

loc:
	cloc --exclude-list-file=.cloc.exclude nuclmm/*.py
	cloc nuclmm/tests/test_*.py

clean:
	rm -rf __pycache__/
