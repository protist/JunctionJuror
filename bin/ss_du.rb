#!/usr/bin/env ruby
#
# Copyright 2013 Lee M. Yeoh (email: "plast-dot-id-dot-au" follows "github")
# This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
# This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# 

require 'optparse'
$options = {}
OptionParser.new do |opts|
  opts.banner='DESCRIPTION:'\
              'TODO
Usage:' # TODO: write description
  opts.on_tail('-h', '--help', 'Show this message') do
    puts opts; exit
  end
  # Verbosity off (or == 0) is base-level information. Verbosity == 1 is
  #   chromosome-level. Verbosity == 2 is junction/gene-level.
  opts.on('-v', '--[no-]verbose [OPT]', 'Run verbosely, optionally at level 2') do |v|
    $options[:verbosity] = (v || 1).to_i
  end
  opts.on('-j', '--junction JUNCTION.BED.LIST',
          'JUNCTION.BED.LIST is a required space-delimited list, with each line '\
          'comprising a path to a junction.bed file, followed by a condition.') do |j|
    $options[:junction_list] = File.expand_path(j)
  end
  opts.on('-g', '--ref_gff GFF_FILE',
          'Reference GFF_FILE (e.g. from eupathdb) is required.') do |g|
    $options[:refgff_path] = File.expand_path(g)
  end
  opts.on('-o', '--output OUTPUT_PATH', 'OUTPUT_PATH is required.') do |o|
    $options[:output_path] = File.expand_path(o)
  end
  opts.on('-t', '--threshold NUMBER', 'NUMBER replicates are needed to '\
      'confirm a junction. Defaults to 2.') do |t|
    $options[:min_replicates] = t.to_i
  end
end.parse!

# Require Ruby 1.9.2 for ordered hash. This is quicker than continual hash.sort.
min_release = '1.9.2'
ruby_release = RUBY_VERSION
if ruby_release < min_release
  abort("This script requires Ruby version #{min_release} or later. You are running #{ruby_release}.")
end

# Mandatory "options".
raise OptionParser::MissingArgument if $options[:junction_list].nil? ||
    $options[:refgff_path].nil? ||
    $options[:output_path].nil?
$options[:verbosity] = 0 if !$options[:verbosity]
$options[:min_replicates] = 2 if !$options[:min_replicates]

# Check options.
to_abort=false
if !(File.file? $options[:junction_list])
  $stderr.puts "Error: junction.bed list #{$options[:junction_list]} doesn't exist."
  to_abort=true
end
if !(File.exists? $options[:refgff_path])
  $stderr.puts "Error: reference gff #{$options[:refgff_path]} doesn't exist."
  to_abort=true
end
if File.exists? $options[:output_path]
  $stderr.puts "Error: output path #{$options[:output_path]} already exists."
  to_abort=true
end
raise OptionParser::InvalidArgument if to_abort

################################################################################
### Define Junction class
class UserJunctions
  # A UserJunctions object stores information from TopHat outputs, specifically
  # the chromosome and the coordinates of the skipped bases in the junction.
  # It will also store the associated gene ID.
  # @junctions
  # {:condition => {:chromosome => [[start, stop, count],…]}}}
  # later,
  # {:condition => {:chromosome => [[start, stop, count, :gene_id],…]}}
  def initialize(junctions = {})
    @junctions = junctions
  end

  attr_reader(:junctions)

  # Create a new condition and chromosome if necessary, then add to count for
  #   coordinates.
  def write_junctions_for_replicate(condition, chromosome, start, stop)
    # Create new condition hash if it doesn't exist.
    @junctions[condition] ||= {}
    # Create new chromosome array if it doesn't exist.
    @junctions[condition][chromosome.to_sym] ||= []
    # Append coordinates to this array if necessary; add count.
    # TODO: This reads the entire array each time. If this is too slow, it might
    #   be better to read each replicate into individual objects, sort, then
    #   collate sequentially.
    junction_index = @junctions[condition][chromosome.to_sym].index do |junction|
      junction[0] == start && junction[1] == stop
    end
    if junction_index
      @junctions[condition][chromosome.to_sym][junction_index][2] += 1
    else
      @junctions[condition][chromosome.to_sym].push [start,stop,1]
    end
  end

  # Return a count of the total number of junctions with specified replicates.
  #   If `replicates` is false, then returns the total count.
  def count_junctions(replicates)
    total = 0
    @junctions.each_value do |junctions_by_chromosome|
      junctions_by_chromosome.each_value do |junctions|
        if replicates
          total += junctions.count {|junction| junction[2] == replicates}
        else
          total += junctions.count
        end
      end
    end
    total
  end

  # Remove junctions according to number of confirming replicates.
  def prune!
    @junctions.each do |condition, junctions_by_chromosome|
      junctions_by_chromosome.each do |chromosome, _|
        @junctions[condition][chromosome].reject! do |junction|
          junction[2] < $options[:min_replicates]
        end
      end
    end
  end

  # Sort junctions according to start coords.
  def sort!
    @junctions.each do |condition, junctions_by_chromosome|
      junctions_by_chromosome.each do |chromosome, _|
        @junctions[condition][chromosome].sort!
      end
    end
  end
end

################################################################################
### Define ReferenceGFF class
class ReferenceGFF
  # A ReferenceGFF object stores gene models, specifically chromosome, gene ID,
  # and terminal start and stop coordinates (e.g. of the CDS).
  def initialize(genes_by_chromosome = {})
    @genes_by_chromosome = genes_by_chromosome
  end

  attr_reader(:genes_by_chromosome) # for debugging

  # Create a new chromosome and gene id if necessary, and replace start and stop
  # coordinates if they increase the boundaries.
  # {:chromosome => {:geneID => [start, stop]} }
  def write_gene(chromosome, gene_id, start, stop)
    # Create new chromosome hash if it doesn't exist.
    @genes_by_chromosome[chromosome.to_sym] ||= {}
    # Create new gene_id hash if it doesn't exist.
    @genes_by_chromosome[chromosome.to_sym][gene_id.to_sym] ||= [start, stop]
    # Compare coordinates to existing ones (redundancy if you've just created
    #   the gene, but I'm not sure if it's still quicker to use ||= )
    if start < @genes_by_chromosome[chromosome.to_sym][gene_id.to_sym].first
      @genes_by_chromosome[chromosome.to_sym][gene_id.to_sym][0] = start
    end
    if stop > @genes_by_chromosome[chromosome.to_sym][gene_id.to_sym].last
      @genes_by_chromosome[chromosome.to_sym][gene_id.to_sym][1] = stop
    end
  end

  # For each chromosome, sort genes by start coordinates.
  def sort!
    @genes_by_chromosome.each do |chromosome, genes|
      ordered_genes = Hash[genes.sort_by { |_, coords| coords }]
      @genes_by_chromosome[chromosome] = ordered_genes
      # I'd prefer the following, but Ruby doesn't have sort_by! for hashes. :(
      #   @genes_by_chromosome[chromosome].sort_by! { |_, coords| coords[0] }
    end
  end

  # Check to see if adjacent genes overlap (or touch). If so, then notify the
  #   user and abort. Users should then manually fix it (e.g. by deleting one of
  #   the genes).
  # N.B. for ToxoDB files, there are only 3 and 14 overlaps for GT1 and ME49
  #   respectively anyway.
  def check_overlaps
    overlap_count = 0
    total_count = 0
    @genes_by_chromosome.each do |chromosome, genes_for_this_chromosome|
      prev_gene_id = nil
      genes_for_this_chromosome.each do |gene_id, coords|
        if prev_gene_id
          total_count += 1
          if genes_for_this_chromosome[prev_gene_id].last >= coords.first - 1
            overlap_count += 1
            puts "#{Time.new}:   ERROR! On #{chromosome}, genes "\
              "#{prev_gene_id.to_s} and #{gene_id.to_s} overlap by (at least) "\
              "#{genes_for_this_chromosome[prev_gene_id].last - coords.first + 1}"\
              ' bp.'
          end
        end
        prev_gene_id = gene_id
      end
    end
    if overlap_count > 0
      abort("#{Time.new}:   ERROR! #{overlap_count} overlaps reported.")
    end
  end
end


################################################################################
### Read junction.bed file from user's experimental data
# Parse junction.bed list.
# Junction paths can be relative to the list path, or absolute. Ignore comments.
puts "#{Time.new}: Parsing junction.bed list."
junction_list = {} # {:condition => ['path/rep1', '/path/rep2'],…}

File.open($options[:junction_list]).each do |line|
  # These are the parts of each line that we need.
  split_line = line.split
  if split_line != [] && /^[^#]/ =~ split_line[0]
    junction_list[split_line[1].to_sym] ||= []
    junction_list[split_line[1].to_sym].push File.expand_path(split_line[0],\
        File.dirname($options[:junction_list])) # relative to junctions.bed.list
  end
end

puts "#{Time.new}:   #{junction_list.count} conditions: "\
     "#{junction_list.keys.collect {|cond| cond.to_s}.join(', ')}"
puts "#{Time.new}:   with replicates: "\
     "#{junction_list.values.collect {|cond| cond.count}.join(', ')}."

# Parse junction.bed files themselves.
#   http://genome.ucsc.edu/FAQ/FAQformat.html#format1
#   chr feature_start-1 feature_end (name) (depth) (strand) (bold_feature_start-1)
#     (bold_feature_end-1) (rgb) (#exons) blocksizes block_start_relative_to_feature_start
#   For TopHat out, ignore (parenthesised):
#     name, strand, bold_features and rgb are irrelevant; depth might be useful
#     to define cutoffs for "real" splicing, but is ignored here since we have
#     replicates; #exons is always 2
#   Start coordinate == $1 + $10:0 + 1 == feature_start-1 + 1st_blocksize + 1.
#   Stop coordinate == $2 - $10:1 == feature_end - 2nd_blocksize
#     or $1 + $11:1 == feature_start-1 + 2nd_block_start_relative_to_feature_start
#   Also need chr.
puts "#{Time.new}: Importing junction.bed files."
junctions = UserJunctions.new
junction_list.each do |cond, paths|
  paths.each do |path|
    File.open(path).drop(1).each do |line| # skip first line
      split_line = line.split("\t")
      block_sizes = split_line[10].split(',')
      abort("Error: expected two blocks for junction #{split_line[3]} in #{path}.") if block_sizes.count != 2
      start = split_line[1].to_i + block_sizes[0].to_i + 1
      stop = split_line[2].to_i - block_sizes[1].to_i
      junctions.write_junctions_for_replicate(cond, split_line[0], start, stop)
    end
  end
end

# Import gff.
# Since the UTRs from eupathdb are not all defined, best to be consistent and
#   define genes as being from first to last CDS (except for tRNA and rRNA).
# In field 9 -> Parent=rna_TGME49_203135-1 (N.B. all CDS are /-1$/).
# N.B. if strand == "-", arranged in reverse numerical order, but let's not make
#   either assumption here, just in case this file does not come from eupathdb.
# I think these files don't contain any alternatively-spliced genes. At least,
#   (for TGGT1 and TGME49) no entries have {$3 == "mRNA"} and contain "-2".
#   Even if they did, this script only considers the terminal exons of all CDS
#   sharing the same parent (without -1).
# GFF header lines won't have split_line[2].
puts "#{Time.new}: Parsing reference GFF file."
refgff = ReferenceGFF.new
File.open($options[:refgff_path]).each do |line|
  split_line = line.split("\t")
  skip = false
  if split_line[2] == 'CDS'
    match_data = /Parent=rna_(TG[^_]{2,4}_\d*[A-Z]?)-1(;|$)/.match(split_line[8])
  elsif split_line[2] == 'tRNA' || split_line[2] == 'rRNA'
    match_data = /Parent=(TG[^_]{2,4}_\d*)(;|$)/.match(split_line[8])
  else
    skip = true
  end
  if !skip
    if match_data
      refgff.write_gene(split_line[0], match_data[1], split_line[3].to_i, split_line[4].to_i)
    else
      abort('ERROR: the following line does not contain a gene_id in the ' \
          "expected format\n#{line}")
    end
  end
end

# Check that the gff is ordered.
puts "#{Time.new}:   Sorting reference gff."
refgff.sort!
p refgff.genes_by_chromosome
puts "#{Time.new}:   Checking reference gff for overlapping genes."
refgff.check_overlaps

# Prune junctions with too few replicates, and output some statistics.
puts "#{Time.new}: Pruning junctions in <#{$options[:min_replicates]} replicates."
pre_count = junctions.count_junctions(false)
junctions.prune!
post_count = junctions.count_junctions(false)
puts "#{Time.new}:   #{pre_count - post_count} junctions removed."
max_replicates = junction_list.values.collect {|cond| cond.count}.max
($options[:min_replicates]..max_replicates).each do |replicates|
  puts "#{Time.new}:   #{junctions.count_junctions(replicates)} junctions "\
      "confirmed in #{replicates} replicates."
end

# Sort junctions and gff.
puts "#{Time.new}: Sorting junctions."
junctions.sort!

# Assign genes to junctions. A gene is assigned if at least one matching base
#   from both flanks overlaps with the CDS. Otherwise it is discarded.
