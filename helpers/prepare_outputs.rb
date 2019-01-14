require 'bio'
require 'pry'

DIR_PATH = "../database/outputs"

def prepare_file(file_name)
  file = Bio::FlatFile.open("#{DIR_PATH}/#{file_name}")
  upcase_file = File.open("#{DIR_PATH}/upcase/#{file_name}", "w")

  file.entries.each do |fasta_entry|
    sequence = fasta_entry.data.delete("\n").upcase
    chars = sequence.split("").uniq.sort
    unless chars == ["-", "A", "C", "G", "T"] || chars == ["A", "C", "G", "T"]
      puts "#{file_name} is invalid"
      puts "\T Chars in sequence: #{sequence.split("").uniq.sort}"
      File.delete(upcase_file.path)
    else
      upcase_file.puts(">" + fasta_entry.definition)
      upcase_file.puts(sequence)
    end
  end
end

system 'mkdir', '-p', "#{DIR_PATH}/upcase/"

Dir.foreach('#{DIR_PATH}') do |file_name|
  if file_name.split('.').last == 'fasta'
    prepare_file(file_name)
  end
end