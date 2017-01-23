#!/bin/bash

#pandoc --reference-links --no-wrap -t rst "README.md" > "README.rst"
#python2 setup.py --long-description | rst2html.py > build/README.html
python2 setup.py sdist upload
