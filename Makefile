test:
	pytest -v --cov=nuclmm --cov-report term-missing nuclmm/*.py nuclmm/tests/*.py

style:
	pep8 nuclmm/*.py nuclmm/tests/*.py

clean:
	rm -rf __pycache__/
