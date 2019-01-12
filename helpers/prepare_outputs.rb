require 'bio'
require 'pry'

def prepare_file(file_name)
  file = Bio::FlatFile.open("../database/outputs_mafft/#{file_name}")
  upcase_file = File.open("../database/outputs_mafft/upcase/#{file_name}", "w")

  file.entries.each do |fasta_entry|
    sequence = fasta_entry.data.delete("\n").upcase
    chars = sequence.split("").uniq.sort
    unless chars == ["-", "A", "C", "G", "T"] || chars == ["A", "C", "G", "T"]
      puts "#{file_name} is invalid"
      puts "\T Chars in sequence: #{sequence.split("").uniq.sort}"
    end

    upcase_file.puts(fasta_entry.definition)
    upcase_file.puts(sequence)
  end
end

system 'mkdir', '-p', "../database/outputs_mafft/upcase/"

Dir.foreach('../database/outputs_mafft') do |file_name|
  if file_name.split('.').last == 'fasta'
    prepare_file(file_name)
  end
end