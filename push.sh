#!/bin/bash
# git pusher will test the program first and generate documentation
doxygen "./doxygen/Doxyfile mac"
git add "./FragGeneScan1.31/doc/" && git commit -m"new doc (auto generated)"
git push origin master
