cut -d '|' -f 2 ../fasta_files/uniprot_ids_homologues.txt > accessions.txt

input_file="$1"
output_dir="../fasta_files/species_files"

mkdir -p $output_dir

for acc in $(cat accessions.txt); do
    curl -s "https://rest.uniprot.org/uniprotkb/$acc.fasta"
done > ../fasta_files/all_blast_query_sequences.fasta

csplit -s -z -f "temp_" "../fasta_files/all_blast_query_sequences.fasta" '/^>/' '{*}'

for temp_file in temp_*; do
    if [[ -s "$temp_file" ]]; then
        first_line=$(head -n1 $temp_file)
        species=$(echo $first_line | grep "^>" | cut -d"=" -f2 | cut -d" " -f1-2 | sort | uniq)
        safe_species=$(echo "$species" | sed 's/ /_/')

        cat $temp_file >> "${output_dir}/${safe_species}.fasta"

    fi
    rm $temp_file
done

echo "Separation done, files at ${output_dir}"