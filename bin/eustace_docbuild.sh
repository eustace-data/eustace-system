#!/bin/bash
#
#  Build EUSTACE documentation and output to standard location
# Latex documentation not in production, for the moment


# retrieve CODE_PATH string from EUSTACE system configuration
SYSTEM_PATH=`python2.7 -c "import eustaceconfig;print eustaceconfig.SYSTEM_PATH"`
DOCS_PATH=`python2.7 -c "import eustaceconfig;print eustaceconfig.DOCS_PATH"`
#DOCS_PATH='/home/users/fcapponi/code/svn-eustace/system/try_doc'

# make paths to sphinx docs and output
INPUTDIR=$SYSTEM_PATH/data/docs/sphinx

OUTPUTDIR=$DOCS_PATH
TO_SAVE_1=".htaccess"
TO_SAVE_2=".htpasswd"
find ${OUTPUTDIR} -type f -not -name ${TO_SAVE_1} -not -name ${TO_SAVE_2} -delete

LATEXDIR=$DOCS_PATH/latex
# build them
mkdir -p $LATEXDIR
sphinx-build -b html $INPUTDIR $OUTPUTDIR

# changing permissions on target directories: this allows other users to modify the documentation
#7-owner: rwx
#7-group: rwx
#5-others: r-x

#find ${OUTPUTDIR}/* -type d -exec chmod 775 {} \;

#sphinx-build -b latex $INPUTDIR $LATEXDIR
#pushd $LATEXDIR

# To be used temporarily: does not fix all the problems
#sed -i 's/LATEX = latexmk/LATEX = pdflatex/g' Makefile 
#sed -i 's/PDFLATEX = latexmk/PDFLATEX = pdflatex/g' Makefile 
#make 
#make 
#popd
