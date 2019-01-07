require 'bio'
require 'pry'

file = Bio::FlatFile.open("../database/HIV1_REF_2010_genome_DNA.fasta")
fasta_entries = file.entries

cleaned_sequences = []

fasta_entries.each do |fasta_entry|
  cleaned_sequences << {
    sequence: fasta_entry.data.delete(" -").delete("\n"),
    definition: fasta_entry.definition 
  }
end

appropriate_sequences = []

cleaned_sequences.each do |entry|
  nucleotides = entry[:sequence].split("").uniq
  if nucleotides.size == 4 && nucleotides.sort == ["A", "C", "G", "T"]
    appropriate_sequences << entry
  end
end

sequence_pairs = []

appropriate_sequences.each_with_index do |entry1, index|
  appropriate_sequences.drop(index + 1).each do |entry2|
    sequence_pairs << { entry1: entry1, entry2: entry2 }
  end
end

sequence_pairs.each_with_index do |pair, index|
  f = File.open("../database/pairs/p#{index+1}.fasta", "w")
  f.puts(">" + pair[:entry1][:definition])
  f.puts(pair[:entry1][:sequence])
  f.puts(">" + pair[:entry2][:definition])
  f.puts(pair[:entry2][:sequence])
end