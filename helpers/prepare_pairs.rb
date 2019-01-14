require 'bio'
require 'pry'

FILENAMES = ['random_100', 'random_5000']

FILENAMES.each do |file_name|
  file = Bio::FlatFile.open("../database/#{file_name}.fasta")
  fasta_entries = file.entries

  cleaned_sequences = []

  fasta_entries.each do |fasta_entry|
    sequence = fasta_entry.data.delete(" -").delete("\n").upcase
    nucleotides = sequence.split("").uniq
    if nucleotides.size == 4 && nucleotides.sort == ["A", "C", "G", "T"]
      cleaned_sequences << {
        sequence: sequence,
        definition: fasta_entry.definition 
      }
    end
  end

  pair_number = 0
  cleaned_sequences.each_with_index do |entry1, index|
    cleaned_sequences.drop(index + 1).each do |entry2|
      system 'mkdir', '-p', "../database/pairs/#{file_name}/"
      pair_number += 1

      f = File.open("../database/pairs/#{file_name}/p#{pair_number}.fasta", "w")
      f.puts(">" + entry1[:definition])
      f.puts(entry1[:sequence])
      f.puts(">" + entry2[:definition])
      f.puts(entry2[:sequence])
    end
  end
end