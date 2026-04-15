.PHONY: test cli-help

test:
	PYTHONPATH=src python3 -m unittest discover -s tests -v

cli-help:
	PYTHONPATH=src python3 -m dna_toolkit.cli --help

