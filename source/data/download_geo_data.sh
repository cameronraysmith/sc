#!/usr/bin/env bash

######################
#
# description:
#
#   download supplementary or other files from a GEO series/sample accession ID
#
# usage:
#
#   $ ./download_geo_data.sh \
#       -a GSE132771 \
#       -f 'ftp.*RAW.*' \
#       -j '..|.supplementary_files?|..|.url?|select(length>0)'
#
#   it is intended to be able to rerun the script with
#   identical parameters without modifying existing data
#   (i.e. it is intended to have idempotent side effects)
#
# dependencies:
#
#   python:
#
#     ffq ( if failing to produce $ACCESSION.json,
#           check issues with latest version )
#
#   system:
#
#     jq aria2c tree
#
######################

help()
{
    # display help
    printf "\ndownload_geo_data:\n\n"
    printf "\tdownload supplementary or other files from a GEO series/sample accession ID\n\n"
    printf "\tUsage:  %s [-h] [-a ACCESSION] [-f URL_FILTER_PATTERN] [-j JSON_FILTER_PATTERN]\n\n" "$0"
    printf "\tExample:  $ ./download_geo_data.sh \\
       \t\t\t-a GSE132771 \\
       \t\t\t-f 'ftp.*RAW.*' \\
       \t\t\t-j '..|.supplementary_files?|..|.url?|select(length>0)'\n\n"
    echo "options:"
    echo "   -h     print this help"
    echo "   -a     set GEO accession (e.g. GSE132771)"
    echo "   -f     set URL search pattern (e.g. 'ftp.*RAW*' )"
    echo "   -j     set JSON search pattern (e.g. '..|.supplementary_files?|..|.url?|select(length>0)' )"
    echo
}

#-- debugging (comment to reduce stderr output)
#-- https://wiki.bash-hackers.org/scripting/debuggingtips
# export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
# set -o xtrace

# Set default variables

# target GEO accession ID
# examples: GSE132771 GSE133549
ACCESSION=GSE171524

# jq string to filter JSON returned by ffq
# test your pattern on stub data at https://jqplay.org/#
# example: '."geo_samples"[]."supplementary_files"[].url'
JSON_FILTER_PATTERN='..|.supplementary_files?|..|.url?|select(length>0)'

# grep regex pattern for filtering urls
# examples: 'NML.*All' 'ftp.*RAW*'
URL_FILTER_PATTERN='ftp.*RAW*'


# Get the options
while getopts ":ha:f:j:" option; do
    case $option in
        h) # display Help
            help
            exit 1 ;;
        a) # Enter a GEO accession
            ACCESSION=$OPTARG ;;
        f) # Enter a filter pattern for URLs to download
            URL_FILTER_PATTERN=$OPTARG ;;
        j) # Enter a filter pattern for JSON returned by GEO API via ffq
            JSON_FILTER_PATTERN=$OPTARG ;;
        \?) # Invalid option
            echo "Error: Invalid option"
            help
            exit 1 ;;
    esac
done


#-- declare variables
DOWNLOAD_URLS_FILE=download_urls.txt # file to save list of urls to be downloaded
DOWNLOAD_THREADS=4 # number of parallel threads for aria2c to use in downloading
DOWNLOAD_SHARDS=2 # number of fragments to split the file into for download
OUTPUT_FILE_DIRECTORY=supplementary # folder for saving downloaded files

#-- create/cd to directory for accession
mkdir -p "$ACCESSION"
cd "$ACCESSION" || exit 1
#-- NOTE: the remainder is executed inside the $ACCESSION subdirectory
#-- and could be contained within a subshell

#-- get accession metadata using ffq
printf "\n* searching for $ACCESSION \n\n"
if [ ! -f "$ACCESSION".json ]; then
    ffq -l 1 -o "$ACCESSION".json "$ACCESSION"
fi

# -- parse ffq output to download supplementary files
# jq -r '."geo_samples"[]."supplementary_files"[].url' $ACCESSION.json | \#
# jq -r '..|.supplementary_files?|..|.url?|select(length>0)' $ACCESSION.json | \#
printf "\n* filtering $ACCESSION metadata\n\n\t-JSON: $JSON_FILTER_PATTERN\n\t-URL: $URL_FILTER_PATTERN\n\n"
jq -r "$JSON_FILTER_PATTERN" "$ACCESSION".json | \
    grep "$URL_FILTER_PATTERN" > $DOWNLOAD_URLS_FILE
printf "\n* files to be downloaded: \n\n"
cat -b $DOWNLOAD_URLS_FILE

#--  download target files from parsed urls in parallel with aria2c and unzip
printf "\n* downloading files\n"
mkdir -p $OUTPUT_FILE_DIRECTORY
aria2c --input-file $DOWNLOAD_URLS_FILE \
    --continue true \
    -x $DOWNLOAD_SHARDS \
    -j $DOWNLOAD_THREADS \
    -d $OUTPUT_FILE_DIRECTORY

#-- untar downloaded tarballs and list zipped files
printf "\n* untar downloaded tarballs\n\n"
# enter subshell for $OUTPUT_FILE_DIRECTORY
(
    cd $OUTPUT_FILE_DIRECTORY || exit 1
    for f in *.tar; do echo "$f" && tar -C ./ -kxf "$f"; done
    printf "\n* remaining compressed files\n\n"
    for f in *.*z; do echo "  $f"; done
    printf "\n  can be decompressed with:\n\n\t \$ yes n | gunzip -k $ACCESSION/$OUTPUT_FILE_DIRECTORY/*.gz\n\n  or similar for xz, etc\n"
)
# yes n | gunzip -k $OUTPUT_FILE_DIRECTORY/*.gz

#-- print directory tree with file sizes
printf "\n* path to and contents of accession directory\n\n"
pwd
tree --du -h .

printf "\n* scroll up to copy command for decompressing downloaded files \n\n"
exit 0
