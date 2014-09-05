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


end

################################################################################
### Read junction.bed file from user's experimental data
# Parse junction.bed list.
# Junction paths can be relative to the list path, or absolute. Ignore comments.
puts "#{Time.new}: parsing junction.bed list."
junction_list = {} # {:condition => ['path/rep1', '/path/rep2'],…}

File.open($options[:junction_list]).each do |line|
  # These are the parts of each line that we need.
  split_line = line.split
  if split_line != [] && /^[^#]/ =~ split_line[0]
    junction_list[split_line[1].to_sym] ||= []
    junction_list[split_line[1].to_sym].push File.expand_path(split_line[0],\
        File.dirname($options[:junction_list])) # relative to junction.bed.list
  end
end

# Parse junction.bed files.
#   chr feature_start-1 feature_end (name) (depth) (strand) (bold_feature_start-1)
#     (bold_feature_end-1) (rgb) (#exons) blocksizes block_start_relative_to_feature_start
#   For TopHat out, ignore (parenthesised):
#     name, strand, bold_features and rgb are irrelevant; depth might be useful
#     to define cutoffs for "real" splicing, but is ignored here since we have
#     replicates; #exons is always 2
#   Start coordinate == $2 + $11:1 + 1 == feature_start-1 + 1st_blocksize + 1.
#   Stop coordinate == $3 - $11:2 == feature_end - 2nd_blocksize
#     or $2 + $12:2 == feature_start-1 + 2nd_block_start_relative_to_feature_start
#   Also need chr.
junctions = UserJunctions.new

# N.B. junction.bed is not fully sorted.

# Import gff.

# Sort junctions and gff.
puts "#{Time.new}: sorting junctions."
#junctions.sort!
