if [ ! -d "../data_archive" ]
then
    echo "Data archive folder does not exist. Creating folder \"data_archive\""
    mkdir ../data_archive
fi
echo "Moving previous data file to data_archive"
now="$(date +'%d_%m_%Y')"
mv ../data/IEDB_data.tsv ../data_archive/"$now"_IEDB_data.tsv
echo "File has been moved and renamed to "$now"_IEDB_data.tsv"
curl "https://downloads.iedb.org/misc/TCRMatch/IEDB_data.tsv" -o ../data/IEDB_data.tsv