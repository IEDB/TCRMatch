#!/bin/bash
# Download the latest TCRMatch database from the IEDB into data/, moving any
# previous copy to data_archive/ as <day>_<month>_<year>_<file>.
#
#   ./update.sh            data/IEDB_data.tsv          (legacy format, the default -d database)
#   ./update.sh --schema2  data/IEDB_data.schema2.tsv  (one row per receptor group, adds
#                                                       receptor_species; needs TCRMatch >= 1.4;
#                                                       use with -d data/IEDB_data.schema2.tsv)
set -euo pipefail
cd "$(dirname "$0")"

file="IEDB_data.tsv"
case "${1:-}" in
    "") ;;
    --schema2) file="IEDB_data.schema2.tsv" ;;
    *) echo "Usage: $0 [--schema2]" >&2; exit 1 ;;
esac
url="https://downloads.iedb.org/misc/TCRMatch/$file"
tmp="../data/.$file.download"

# Download first, so a failed download never leaves you without a database.
echo "Downloading $url"
if ! curl -fsSL "$url" -o "$tmp"; then
    rm -f "$tmp"
    echo "Download failed; data/$file was not changed." >&2
    exit 1
fi
if [ "$(head -c 11 "$tmp")" != "trimmed_seq" ]; then
    rm -f "$tmp"
    echo "Downloaded file is not a TCRMatch database; data/$file was not changed." >&2
    exit 1
fi

if [ -f "../data/$file" ]; then
    if [ ! -d "../data_archive" ]; then
        echo "Data archive folder does not exist. Creating folder \"data_archive\""
        mkdir ../data_archive
    fi
    now="$(date +'%d_%m_%Y')"
    mv "../data/$file" "../data_archive/${now}_$file"
    echo "Moved previous data/$file to data_archive/${now}_$file"
fi
mv "$tmp" "../data/$file"
echo "Updated data/$file"
