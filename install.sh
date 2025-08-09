#!/bin/bash

#--------------------------------------------------------------------
# a _suite's_ 'install.sh' script installs the MDI and any external suite
# dependencies when running in suite-centric suite mode
# it may install R Shiny packages, but never installs conda environments
#--------------------------------------------------------------------
MDI_CENTRIC="mdi-centric"
SUITE_CENTRIC="suite-centric"
export SUPPRESS_MDI_BASHRC=TRUE
if [[ "$1" = "--yes" || "$1" = "-y" ]]; then
    MDI_FORCE_GIT=TRUE
    MDI_FORCE_APPS=TRUE
fi

#----------------------------------------------------------------------
# base directory of this tool suite, after user cloned it from GitHub
#----------------------------------------------------------------------
export SUITE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export SUITE_NAME=`basename $SUITE_DIR`

#----------------------------------------------------------------------
# discover whether this suite is in an mdi-centric or suite-centric installation
# if it is contained in a parent mdi/suites/* folder, then it is always mdi-centric
# otherwise, the mdi will be installed and run from within this suite directory (i.e., suite-centric)
#----------------------------------------------------------------------
export SUITE_MODE=$MDI_CENTRIC
export MDI_DIR=$SUITE_DIR/../../..
export MDI_TARGET=$MDI_DIR/mdi
if [[ ! -d  "$MDI_DIR/frameworks" || ! -d  "$MDI_DIR/suites" || ! -f  "$MDI_TARGET" ]]; then
    export SUITE_MODE="$SUITE_CENTRIC" # i.e., this is the top-level of a single-suite installation
    export MDI_DIR="$SUITE_DIR/mdi"
    export MDI_TARGET="$MDI_DIR/mdi"
fi
IS_READLINK=`which readlink 2>/dev/null`
if [ "$IS_READLINK" != "" ]; then export MDI_DIR=`readlink -f $MDI_DIR` ; fi

#----------------------------------------------------------------------
# if mdi-centric, exit quietly with a helpful message and no action taken
#----------------------------------------------------------------------
if [ "$SUITE_MODE" = "$MDI_CENTRIC" ]; then
    echo -e "\nNothing to do.\n"
    echo -e "This copy of '$SUITE_NAME' is part of an MDI installation in directory:\n    $MDI_DIR\n"
    echo -e "Use it by calling the 'mdi' command line utility in that directory (multi-suite installation)"
    echo -e "or the 'run' script in the parent of that directory (single suite installation).\n"
    exit 1
fi

#----------------------------------------------------------------------
# get permission to install in suite-centric mode
#----------------------------------------------------------------------
echo
echo "----------------------------------------------------------------------"
echo "Welcome to the Michigan Data Interface (MDI) installer"
echo "----------------------------------------------------------------------"
echo
echo "This script will first initialize '$SUITE_NAME' Stage 1 Pipelines for use"
echo "by cloning the MDI manager utility plus any additional suite dependencies"
echo "(if needed, we will prompt again later regarding Stage 2 Apps)."
echo
if [ "$MDI_FORCE_GIT" != "" ]; then
    PERMISSION=y
else 
    read -p "Do you wish to continue? (type 'y' for 'yes'): " PERMISSION
fi
if [ "$PERMISSION" != "y" ]; then exit; fi

#----------------------------------------------------------------------
# install/update the MDI nested within this suite's directory
#----------------------------------------------------------------------
if [ -d $MDI_DIR ]; then
    cd $MDI_DIR
    git checkout main
    git pull
else 
    git clone https://github.com/MiDataInt/mdi.git
    cd $MDI_DIR
fi
$MDI_DIR/install.sh 1

#----------------------------------------------------------------------
# when installing into a container, copy over any existing base library
#----------------------------------------------------------------------
if [ -e $SUITE_DIR/mdi-tmp-library ]; then
    rm -rf $MDI_DIR/library/*
    mv $SUITE_DIR/mdi-tmp-library/* $MDI_DIR/library
    rm -rf $SUITE_DIR/mdi-tmp-library
fi

#----------------------------------------------------------------------
# programmatically write the suite-centric MDI installation's 'config/suites.yml' file
#----------------------------------------------------------------------
echo "writing mdi/config/suites.yml"
cat \
<(echo -e 'suites:') \
<(
    grep -A 10 -P '^\[remote "origin"\]' $SUITE_DIR/.git/config | 
    grep -v -P '^\[remote "origin"\]' | 
    grep -B 100 -P '^\[' --max-count 1 |
    grep 'url' | 
    sed 's/\s*url\s*=\s*//' | 
    awk '{print "  - "$0}'
) \
<(
    grep -v -P '^\s*#' $SUITE_DIR/_config.yml | 
    grep -P '\S' |
    grep -A 100 -P "^suite_dependencies:" --max-count 1 | 
    grep -v -P "^suite_dependencies:" | 
    grep -B 100 -P '^\S' --max-count 1 | 
    grep -P '\s+-\s+' | 
    sed 's/\s*-\s*//' | 
    awk '{print "  - "$0}'
) > $MDI_DIR/config/suites.yml

#----------------------------------------------------------------------
# re-run the mdi installation to add all required tool suites (including this one)
# installer checks out frameworks to 'latest' and this suite to:
#   latest, if called by user at command line
#   SUITE_VERSION, if called during a container build
#----------------------------------------------------------------------
$MDI_DIR/install.sh 1

#----------------------------------------------------------------------
# clean up this suite root directory by removing most git repo files
# as they are now re-cloned within the nested MDI installation
# leave the utilities to re-install and run the tool suite
#----------------------------------------------------------------------
echo "cleaning up this suite directory"
cd $SUITE_DIR
rm -fr docs .git pipelines shared shiny templates
rm -f  aws-mdi.md .gitignore index.html .lintr overview.md

# -----------------------------------------------------------------------
# function to check for valid singularity
# -----------------------------------------------------------------------
function set_singularity_version_ {
    export SINGULARITY_VERSION=`singularity --version 2>/dev/null | grep -P '^singularity.+version.+'`
}
function set_singularity_version {
    set_singularity_version_ # use system singularity if present
    if [ "$SINGULARITY_VERSION" = "" ]; then # otherwise attempt to load it
        CONFIG_FILE=$MDI_DIR/config/singularity.yml
        if [ -f $CONFIG_FILE ]; then
            LOAD_COMMAND=`grep -P '^load-command:\s+' $CONFIG_FILE | sed -e 's/\"//g' -e 's/load-command:\s*//' | grep -v null | grep -v '~'`
            if [ "$LOAD_COMMAND" != "" ]; then
                $LOAD_COMMAND > /dev/null 2>&1
                set_singularity_version_
            fi
        fi
    fi 
}

#----------------------------------------------------------------------
# discover whether any of the tool suites offer Stage 2 apps
# if so, and not supported by suite-level container, offer to install the apps server 
#----------------------------------------------------------------------
NESTED_SUITES_DIR=$MDI_DIR/suites/definitive
IS_APPS=`perl -e '
foreach my $dir(glob("'$NESTED_SUITES_DIR'/*/shiny/apps/*")){
    use File::Basename;
    -d $dir or next;
    my $app = basename($dir);
    $app =~ m/^_/ and next;
    print "TRUE\n";
    exit;
}'`
if [ "$SUITE_NAME" = "mdi-singularity-base" ]; then IS_BASE_INSTALL=true; fi
if [[ "$IS_APPS" = "TRUE" || "$IS_BASE_INSTALL" != "" ]]; then 
    echo -e "\n----------------------------------------------------------------------"
    echo "'$SUITE_NAME' additionally offers interactive Stage 2 Apps"
    echo "----------------------------------------------------------------------"

    # determine if both system and a suite-level container support apps, no need for direct install
    set_singularity_version
    SUPPORTS_CONTAINERS=`grep -A10 -P '^container:' $SUITE_DIR/_config.yml | grep -P '^\s+supported:\s+true'`
    CONTAINER_HAS_APPS=`grep -A10 -P '^container:' $SUITE_DIR/_config.yml | grep -P '^\s+apps:\s+true'`
    if [[ "$SINGULARITY_VERSION" != "" && "$SUPPORTS_CONTAINERS" != "" && "$CONTAINER_HAS_APPS" != "" && "$IS_BASE_INSTALL" = "" ]]; then
        echo -e "\n'$SUITE_NAME' apps server runs in a suite-level Singularity container."
        echo -e "Further installation not required."

    # if not, proceed with direct apps server installation if confirmed
    else
        echo -e "\nNOTE: The apps server installation requires R and takes many minutes to complete."
        echo -e "Please decline if you will only use Stage 1 data analysis pipelines from this installation.\n"
        if [ "$MDI_SKIP_APPS" != "" ]; then
            PERMISSION=n
        elif [ "$MDI_FORCE_APPS" != "" ]; then
            PERMISSION=y
        else 
            read -p "Install the R Shiny Stage 2 Apps server? (type 'y' for 'yes'): " PERMISSION
        fi
        if [ "$PERMISSION" = "y" ]; then 
            $MDI_DIR/install.sh 2
        else 
            echo "Skipping Stage 2 Apps installation."
        fi
    fi
else 
    echo -e "\n'$SUITE_NAME' does not offer Stage 2 Apps."
fi

# finish up all successful installation paths
echo -e "\nInstallation complete.\n"
