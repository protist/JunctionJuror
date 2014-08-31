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
    $options[:junction_list] = j
  end
  opts.on('-g', '--ref_gff GFF_FILE',
          'Reference GFF_FILE (e.g. from eupathdb) is required.') do |g|
    $options[:refgff_path] = g
  end
  opts.on('-o', '--output OUTPUT_PATH',
          'OUTPUT_PATH is required.') do |o|
    $options[:output_path] = o
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
  puts "Error: junction.bed list #{$options[:junction_list]} doesn't exist."
  to_abort=true
end
if !(File.exists? $options[:refgff_path])
  puts "Error: reference gff #{$options[:refgff_path]} doesn't exist."
  to_abort=true
end
if File.exists? $options[:output_path]
  puts "Error: output path #{$options[:output_path]} already exists."
  to_abort=true
end
if to_abort
  abort
end

# Junction path can be absolute or relative to the list.
