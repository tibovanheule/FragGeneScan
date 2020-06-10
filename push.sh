#!/bin/bash

# git pusher will test the program first and generate documentation
git pull
# generate doxygen
doxygen "./doxygen/Doxyfile mac"
# add new doc files to a commmit
git add "./src/doc/" && git commit -m"new doc (auto generated)"
# navigate to latex doc folder
cd ./src/doc/latex/
# generate new pdf
pdflatex -halt-on-error -output-directory=../../../ -jobname="Fraggenescan documentation" refman.tex
cd ../../../
# commit new pdf
git add "Fraggenescan documentation.pdf" && git commit -m"new doc pdf (auto generated)"
# push all
git push origin master
