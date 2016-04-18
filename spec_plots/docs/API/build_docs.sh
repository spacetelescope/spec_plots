#!/bin/sh
# This works fine.
 sphinx-build -T -E -b html -d _build/doctrees-readthedocs -D language=en . _build/html

# This is what RTD does, need to get this working on my local machine.
#sphinx-build -T -E -b doctest -d _build/doctrees-readthedocs -D language=en . _build/html