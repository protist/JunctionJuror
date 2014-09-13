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
# This script attempted to identify alternative splicing based on junction
#   reads. It doesn't acknowledge intron retention.
#
# Firstly, junction beds files are read from a list, discarding junctions that
#   are confirmed in less than NUMBER replicates. It then assigns junctions to
#   genes as per GFF_FILE, and discards those that do not have at least one
#   matching base within the CDS in both flanks. Then, it detects when multiple
#   junctions overlap, including when adjacent skipped regions abut. Finally,
#   a list of genes with overlapping junctions is output.

require 'optparse'
require 'set'
$options = {}
OptionParser.new do |opts|
  opts.banner='This script reads in junctions.beds files, discarding '\
              'junctions confirmed in less than NUMBER replicates. It assigns'\
              'junctions to genes as per GFF_FILE, discarding intergenic '\
              'junctions. It detects when multiple junctions overlap, and '\
              "outputs the list of genes where this occurs.\nUsage:"
  opts.on_tail('-h', '--help', 'Show this message') do
    puts opts; exit
  end
  # Verbosity off (or == 0) is base-level information. Verbosity == 1 is
  #   chromosome-level. Verbosity == 2 is junction/gene-level.
  $options[:verbosity] = 0
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
  $options[:min_replicates] = 2
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
  # @junctions format:
  #   {:conditionX => {:chrX => [{:coords => [start, stop], :count => num},…]}}
  #   later, add :gene_id => id
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
      junction[:coords].first == start && junction[:coords].last == stop
    end
    if junction_index
      @junctions[condition][chromosome.to_sym][junction_index][:count] += 1
    else
      @junctions[condition][chromosome.to_sym].push({coords: [start,stop], count: 1})
    end
  end

  # Return a count of the total number of junctions with specified replicates.
  #   If `replicates` is zero, then returns the total count.
  def count_junctions(replicates)
    total = 0
    @junctions.each_value do |junctions_by_chromosome|
      junctions_by_chromosome.each_value do |junctions|
        if replicates == 0
          total += junctions.count
        else
          total += junctions.count { |junction| junction[:count] == replicates }
        end
      end
    end
    total
  end

  # Remove junctions according to number of confirming replicates.
  def prune!
    @junctions.each do |condition, junctions_by_chromosome|
      junctions_by_chromosome.each_key do |chromosome|
        @junctions[condition][chromosome].reject! do |junction|
          junction[:count] < $options[:min_replicates]
        end
      end
    end
  end

  # Sort junctions according to start coords.
  #   sort_by is probably more efficient (or similar) than sort here.
  #   https://gist.github.com/protist/380ec7f0a3a53c7835f0
  def sort!
    @junctions.each do |condition, junctions_by_chromosome|
      junctions_by_chromosome.each_key do |chromosome|
        @junctions[condition][chromosome].sort_by! { |junction| junction[:coords] }
      end
    end
  end

  # Remove junctions on chromosomes not found in refgff. Returns array of
  #   missing chromosomes.
  def check_chromosomes!(valid_chromosomes)
    rejects = []
    @junctions.each_value do |junctions_by_chromosome|
      junctions_by_chromosome.reject! do |chromosome_name, _|
        if !(valid_chromosomes.index chromosome_name)
          rejects << chromosome_name
          true
        end
      end
    end
    rejects
  end

  # Assign genes to junctions. A gene is assigned if at least one matching base
  #   from both flanks overlaps with the CDS. Otherwise it is discarded.
  # It's faster to access the gene array with an incrementing counter rather
  #   than conversion to an enumerator, then enum.next (and enum.peek).
  #   https://gist.github.com/protist/f6d344cdc59571947677
  # It's also faster to use map!.compact! rather than inject or flat_map
  #   https://gist.github.com/protist/f153d0465cae3f474081
  # Returns the number of genes that have junctions mapped to them.
  def assign_genes!(refgff)
    matching_gene_list = Set.new
    @junctions.each_value do |junctions_by_chromosome|
      junctions_by_chromosome.each do |chromosome, junctions|
        gene_index = 0
        number_of_genes_to_test = refgff.length(chromosome)
        junctions.map! do |junction|
          parent_gene = nil
          while !parent_gene && (gene_index < number_of_genes_to_test)
            testing_gene = refgff.gene(chromosome, gene_index)
            if (junction[:coords].first - 1 >= testing_gene.last[0]) &&
                (junction[:coords].first - 1 <= testing_gene.last[1]) # Starts in gene.
              if junction[:coords].last + 1 <= testing_gene.last[1] # Stops in gene.
                parent_gene = testing_gene.first
              else # Junction stops post-gene -> break.
                parent_gene = nil
                break
              end
            elsif junction[:coords].first - 1 < testing_gene.last[0] # Starts before gene.
              parent_gene = nil
              break
            else # Junction starts after the gene -> look again.
              gene_index += 1
            end
          end
          if parent_gene.nil?
            nil
          else
            matching_gene_list.add(parent_gene)
            junction.merge({gene_id:parent_gene}) # Add gene_id to junction.
          end
        end.compact! # Remove nils, i.e. those not matching to genes.
      end
    end
    matching_gene_list.count # TODO: This only works for one replicate.
  end

  # Determine when multiple junctions overlap with each other. This is defined
  #   as including abutting junctions, i.e. with skipped bases adjacent.
  #   Return a list of genes. Could also return a list of junctions, but the
  #   specifics of which overlap is even more dependent on read depth. (The only
  #   plausible usage for a list of junctions is comparison between conditions.)
  def genes_with_overlaps
    @as_gene_list = []
    @junctions.each_value do |junctions_by_chromosome|
      junctions_by_chromosome.each_value do |junctions|
        prev_gene_id = nil
        current_junctions_coords = nil
        junctions.each do |junction|
          current_gene_id = junction[:gene_id]
          if current_gene_id != prev_gene_id # New gene.
            current_junctions_coords = [junction[:coords]]
            prev_gene_id = current_gene_id
          elsif current_gene_id != @as_gene_list.last # i.e. not previously added
            current_junctions_coords.each do |checked_junction|
              if (junction[:coords].last >  checked_junction.first - 1) &&
                  (junction[:coords].first < checked_junction.last + 1) # overlap
                @as_gene_list << current_gene_id
                break
              end
            end
            current_junctions_coords << junction[:coords]
          end
        end
      end
    end
    @as_gene_list # TODO: This only works for one replicate.
  end

  # Return the number of genes containing overlapping junctions.
  def count_genes_with_overlaps
    @as_gene_list.count
  end

  # Write genes with overlapping junctions to file.
  def write_genes_to_file(path)
    File.open(path, 'w') do |output_file|
      output_file.puts(@as_gene_list.join("\n"))
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

  # Return list of chromosomes.
  def chromosomes
    @genes_by_chromosome.keys
  end

  # Return number of genes in a given chromosome.
  def length(chromosome)
    @genes_by_chromosome[chromosome].length
  end

  # Return gene for a chromosome and gene-position index in an array.
  def gene(chromosome, gene_index)
    [@genes_by_chromosome[chromosome].keys[gene_index],
     @genes_by_chromosome[chromosome].values[gene_index]]
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
     "#{junction_list.keys.collect { |cond| cond.to_s }.join(', ')}."
puts "#{Time.new}:   with replicates: "\
     "#{junction_list.values.collect { |cond| cond.count }.join(', ')}."

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

# Import GFF.
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

# Sort GFF.
puts "#{Time.new}:   Sorting reference gff."
refgff.sort!
puts "#{Time.new}:   Checking reference gff for overlapping genes."
refgff.check_overlaps

# Prune junctions with too few replicates, and output some statistics.
puts "#{Time.new}: Pruning junctions in < #{$options[:min_replicates]} replicates."
pre_count = junctions.count_junctions(0)
junctions.prune!
post_count = junctions.count_junctions(0)
puts "#{Time.new}:   #{pre_count - post_count} junctions removed."
max_replicates = junction_list.values.collect { |cond| cond.count }.max
total = ($options[:min_replicates]..max_replicates).inject(0) do |sum, replicates|
  junctions_count = junctions.count_junctions(replicates)
  puts "#{Time.new}:   #{junctions_count} junctions "\
      "confirmed in #{replicates} replicates."
  sum + junctions_count
end
puts "#{Time.new}:   #{total} junctions in total."

# Sort junctions and gff.
puts "#{Time.new}: Sorting junctions."
junctions.sort!

# Remove junctions on chromosomes not found in refgff.
puts "#{Time.new}: Verifying junction chromosomes."
rejects = junctions.check_chromosomes!(refgff.chromosomes)
if rejects != []
  puts "#{Time.new}:   WARNING: junction chromosome#{'s' if rejects.count > 1}"\
      ' not found in reference GFF.'
  rejects.each { |reject| puts "#{Time.new}:     #{reject}" }
end

# Assign genes to junctions; discard intergenic junctions. Report on how many
#   junctions are assigned and rejected, and how many genes are associated.
puts "#{Time.new}: Assigning genes to junctions."
pre_count = junctions.count_junctions(0)
matching_gene_count = junctions.assign_genes!(refgff)
post_count = junctions.count_junctions(0)
puts "#{Time.new}:   #{pre_count - post_count} intergenic junctions removed."
puts "#{Time.new}:   #{post_count} junctions remain."
puts "#{Time.new}:   #{matching_gene_count} genes associated with junctions."

# Determine when multiple junctions overlap with each other.
puts "#{Time.new}: Identifying overlapping junctions."
junctions.genes_with_overlaps
puts "#{Time.new}:   #{junctions.count_genes_with_overlaps} genes have "\
    'overlapping junctions.'

# Write genes to file.
puts "#{Time.new}: Writing genes to #{$options[:output_path]}."
junctions.write_genes_to_file($options[:output_path])
