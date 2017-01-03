test:
	py.test -v --cov=nuclmm.py --cov-report term-missing nuclmm.py

pep8:
	pep8 nuclmm.py

clean:
	rm -rf __pycache__/
