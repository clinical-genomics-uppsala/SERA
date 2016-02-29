#/bin/bash

usage() {
cat << EOF
usage: $0	
	
   OPTIONS:
      -h     Show this message
      -p     Project path for output
      -i     Path to inputfile
      -u     Activates graphical user interface (requires Dialog)
      -f     Overwrites all previous output files when encountered

EOF
kill -SIGINT $$
}
#      -c     Path to CNV inputfile (optional)

function fileioerror {
	echo -e "\nWarning:";
	echo -e "The specified file ${1} does not exist.";
}

function settingserror {
	echo -e "\nWarning:";
	echo -e "${1}";
}

# function to check whether a file doesn't exist or should be overwritten
#function checkWrite {
#	if [[ ! -e ${1} || ! -z $FORCE ]]; then
#		return 1;
#	else
#		return 0;
#	fi
#}
