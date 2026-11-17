#!/bin/bash
# Download the latest TCRMatch databases from the IEDB into data/, moving any
# previous copy to data_archive/ as <day>_<month>_<year>_<file>:
#   IEDB_data.schema2.tsv   all IEDB receptors (the default -d database)
#   CEDAR_data.schema2.tsv  receptors from CEDAR (cancer) references
#                           (use with -d data/CEDAR_data.schema2.tsv)
#
# Both have one line per receptor group and a receptor_species column, and
# need TCRMatch 1.4 or later. The legacy IEDB_data.tsv and CEDAR_data.tsv of
# TCRMatch 1.3 and earlier are no longer supported.
set -euo pipefail
cd "$(dirname "$0")"

if [ $# -gt 0 ]; then
    echo "Usage: $0" >&2
    exit 1
fi

update() {
    local file="$1"
    local url="https://downloads.iedb.org/misc/TCRMatch/$file"
    local tmp="../data/.$file.download"

    # Download first, so a failed download never leaves you without a database.
    echo "Downloading $url"
    if ! curl -fsSL "$url" -o "$tmp"; then
        rm -f "$tmp"
        echo "Download failed; data/$file was not changed." >&2
        return 1
    fi
    if [ "$(head -c 11 "$tmp")" != "trimmed_seq" ] || ! head -1 "$tmp" | grep -q "receptor_species"; then
        rm -f "$tmp"
        echo "Downloaded file is not a TCRMatch 1.4 database; data/$file was not changed." >&2
        return 1
    fi

    if [ -f "../data/$file" ]; then
        if [ ! -d "../data_archive" ]; then
            echo "Data archive folder does not exist. Creating folder \"data_archive\""
            mkdir ../data_archive
        fi
        local now
        now="$(date +'%d_%m_%Y')"
        mv "../data/$file" "../data_archive/${now}_$file"
        echo "Moved previous data/$file to data_archive/${now}_$file"
    fi
    mv "$tmp" "../data/$file"
    echo "Updated data/$file"
}

status=0
for file in IEDB_data.schema2.tsv CEDAR_data.schema2.tsv; do
    update "$file" || status=1
done
exit $status
