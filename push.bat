rem git pusher will test the program first and generate documentation

rem generate doxygen
doxygen "./doxygen/Doxyfile mac"
rem add new doc files to a commmit
git add "./src/doc/" && git commit -m"new doc (auto generated)"
rem navigate to latex doc folder
cd ./src/doc/latex/
rem generate new pdf
pdflatex -halt-on-error -output-directory=../../../ -jobname="Fraggenescan documentation" refman.tex
cd ../../../
rem commit new pdf
git add "Fraggenescan documentation.pdf" && git commit -m"new doc pdf (auto generated)"
rem push all
git push origin master
pause
