#!/usr/bin/env ruby

# == Synopsis
#
# Convert a sequence of a given alphabet to ASCII characters corresponding
# to the alphabet character numbers.
#
# == Usage
#
# convert_seq.rb [ -h | --help ] [ -H | --discard-header ] [ -s | --case-sensitive] [alphabet] [ filename ]
#
# The alphabet is given as a string of all characters in the order they should
# be represented as upon conversion.  All instances of the first character will
# be replaced with ASCII character 0, the second with 1, and so on.  More than
# 256 characters in an alphabet are not supported.  Case insensitivity is
# assumed unless the appropriate flag is given to override this behavior.
# If an alphabet is not specified, "ACGT" is assumed.
#
# == Author Todd Wasson

require 'rdoc/ri/ri_paths'
require 'rdoc/usage'
require 'optparse'

def parse_command_line
  args_hash = Hash.new
  opts = OptionParser.new
  opts.on("--help", "-h") {RDoc::usage}
  opts.on("--discard-header", "-H") {args_hash['discard-header'] = true}
  opts.on("--case-sensitive", "-s") {args_hash['case-sensitive'] = true}
  opts.on("--alphabet", "-a", "=STRING", String) {|val| args_hash['alphabet'] = val}

  argv_rest = opts.parse(ARGV) rescue RDoc::usage("usage")
  argv_rest = nil if argv_rest.length == 0
  if !(args_hash['filename'] = argv_rest.shift)
    args_hash['stdin'] = true
  end

  if !args_hash['alphabet']
    args_hash['alphabet'] = "ACGT"
  end

  [argv_rest, args_hash]
end

argv_rest, args_hash = parse_command_line

alphabet_map = Hash.new
count = 0
args_hash['alphabet'].split(//).each {|x| alphabet_map[x[0]] = count.chr; count += 1}

f = args_hash['filename'] ? File.open(args_hash['filename']) : STDIN
d = f.readlines

d.shift if args_hash['discard-header']
d = d.collect{|x| x.chomp.strip}.join
d.upcase! unless args_hash['case-sensitive']
d.each_byte {|b| print alphabet_map[b]}
