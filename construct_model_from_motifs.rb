#!/usr/bin/env ruby

# == Synopsis
#
# Construct a model for use with binding model including background (unbound)
# sequence and optionally nucleosome and a number of TF motifs.
#
# == Usage
#
# construct_hmm_from_motifs.rb [ dbf-ids ] [ -n | --include-nucleosome ]
# [ -h | --help ]
#
# DBF IDs should be space-delimited.
#
# == Author Todd Wasson

require 'construct_nucleosome_model'
require 'macisaac_motif_file_utils'
require 'pbm_motif_file_utils'
require 'rdoc/ri/ri_paths'
require 'rdoc/usage'
require 'optparse'

def parse_command_line
  args_hash = Hash.new
  opts = OptionParser.new
  opts.on("--help", "-h") {RDoc::usage}
  opts.on("--include-nucleosome", "-n") {args_hash[:include_nucleosome] = true}

  argv_rest = opts.parse(ARGV) rescue RDoc::usage("usage")
  args_hash[:gene_ids] = argv_rest.collect{|x| x.upcase}
  argv_rest = nil if argv_rest.length == 0

  if (!args_hash[:include_nucleosome] && (args_hash[:gene_ids].length == 0))
    RDoc::usage("usage")
    exit
  end

  [argv_rest, args_hash]
end




def parse_harbison_motif_file(filename)
  h = Hash.new
  d = IO.readlines(filename).join

  # specially crafted magical regex to parse the specific format of this file
  d.scan(/Log-odds.*?Probability matrix.*?\n.*?\n(.*?)\s*\n(.*?)\s*\n(.*?)\s*\n(.*?)\s*\n.*?Source:\s+(.*?)\s/m).each {|i|
    rows = i[0..3].collect{|j| foo = j.split(/\s+/); foo.shift; foo}
    gene = i[-1]
    h[gene] = []
    rows[0].length.times { h[gene].push([]) }
    rows.each {|j|
      j.each_index {|k|
	h[gene][k].push(j[k].to_f)
      }
    }
  }

  h
end

def build_TATA_motif(motifs_alphabet_map)
#  seqs = ["TATATAAA", "TATATAAG", "TATAAAAA", "TATAAAAG", "TATAAATA", "TATATATA"]
#  weights = [1, 1, 1, 1, 1, 1]
  seqs = ["TATAAAAA", "TATAAAAG", "TATAAATA", "TATAAATG", "TATATAAA", "TATATAAG", "TATATATA", "TATATATG"]
  weights = [3168, 2096, 2720, 984, 4400, 2160, 6016, 2320]
  sum = weights.inject( nil ) { |sum,x| sum ? sum+x : x }
  freqs = []
  seqs[0].length.times { freqs.push([0, 0, 0, 0])}

  for i in 0..(seqs[0].length-1)
    for j in 0..(seqs.length-1)
      freqs[i][motifs_alphabet_map[seqs[j][i].chr]] += weights[j].to_f/sum
    end
  end
  freqs
end

def build_ORC_motif(motifs_alphabet_map, background_motif, padding = 0)
  bg = [0, 0, 0, 0]
  bg[motifs_alphabet_map["A"]] = background_motif[0]
  bg[motifs_alphabet_map["C"]] = background_motif[1]
  bg[motifs_alphabet_map["G"]] = background_motif[2]
  bg[motifs_alphabet_map["T"]] = background_motif[3]
# Aparicio motif
  a = [0.380952380952381, 0.333333333333333, 0.309523809523810, 0.380952380952381, 0.0238095238095238, 0.142857142857143, 0.0714285714285714, 0.904761904761905, 0.0714285714285714, 0.5, 0.0, 0.0, 0.0476190476190476, 0.452380952380952, 0.0, 0.0714285714285714, 0.190476190476190, 0.357142857142857, 0.166666666666667, 0.380952380952381, 0.452380952380952, 0.476190476190476, 0.523809523809524, 0.404761904761905, 0.309523809523810, 0.285714285714286, 0.190476190476190, 0.428571428571429, 0.380952380952381, 0.404761904761905, 0.285714285714286, 0.190476190476190, 0.285714285714286]
  c = [0.0714285714285714, 0.0714285714285714, 0.0238095238095238, 0.0, 0.0, 0.0476190476190476, 0.0476190476190476, 0.0238095238095238, 0.214285714285714, 0.0238095238095238, 0.0, 0.0476190476190476, 0.0, 0.0, 0.119047619047619, 0.0714285714285714, 0.0476190476190476, 0.190476190476190, 0.238095238095238, 0.119047619047619, 0.0952380952380952, 0.142857142857143, 0.0714285714285714, 0.0952380952380952, 0.166666666666667, 0.238095238095238, 0.261904761904762, 0.166666666666667, 0.214285714285714, 0.190476190476190, 0.0714285714285714, 0.0, 0.119047619047619]
  g = [0.0476190476190476, 0.0714285714285714, 0.0238095238095238, 0.0476190476190476, 0.0238095238095238, 0.0, 0.0, 0.0238095238095238, 0.0, 0.404761904761905, 0.0, 0.0, 0.0, 0.0, 0.714285714285714, 0.214285714285714, 0.0952380952380952, 0.190476190476190, 0.166666666666667, 0.214285714285714, 0.166666666666667, 0.119047619047619, 0.166666666666667, 0.0952380952380952, 0.214285714285714, 0.166666666666667, 0.214285714285714, 0.142857142857143, 0.261904761904762, 0.0714285714285714, 0.0714285714285714, 0.0238095238095238, 0.0714285714285714]
  t = [0.5, 0.523809523809524, 0.642857142857143, 0.571428571428571, 0.952380952380952, 0.80952380952381, 0.880952380952381, 0.0476190476190476, 0.714285714285714, 0.0714285714285714, 1.0, 0.952380952380952, 0.952380952380952, 0.547619047619048, 0.166666666666667, 0.642857142857143, 0.666666666666667, 0.261904761904762, 0.428571428571429, 0.285714285714286, 0.285714285714286, 0.261904761904762, 0.238095238095238, 0.404761904761905, 0.309523809523810, 0.309523809523810, 0.333333333333333, 0.261904761904762, 0.142857142857143, 0.333333333333333, 0.571428571428571, 0.785714285714286, 0.523809523809524]

  freqs = []
  a.length.times { freqs.push([0, 0, 0, 0])}
  for i in 0..(a.length-1)
    freqs[i][motifs_alphabet_map["A"]] = a[i]
    freqs[i][motifs_alphabet_map["C"]] = c[i]
    freqs[i][motifs_alphabet_map["G"]] = g[i]
    freqs[i][motifs_alphabet_map["T"]] = t[i]
  end
  padding.times{freqs.unshift(bg)}
  padding.times{freqs.push(bg)}
  freqs
end


argv_rest, args_hash = parse_command_line

reverse_map = {"A" => "T", "T" => "A", "C" => "G", "G" => "C"}
model_alphabet = ["A", "C", "G", "T"]
background_motif = [0.308512, 0.191488, 0.191488, 0.308512]

#motifs_alphabet_map = {"A" => 0, "C" => 1, "T" => 2, "G" => 3} # the tamo file has these in a weird order
#motifs = parse_macisaac_motif_file("macisaac_2006_motifs.tamo")

motifs_alphabet_map = {"A" => 0, "C" => 1, "G" => 2, "T" => 3}
motifs = parse_pbm_motif_directory("pbm")

motifs["TBP"] = build_TATA_motif(motifs_alphabet_map)
motifs["ORC"] = build_ORC_motif(motifs_alphabet_map, background_motif)

motifs_to_include = args_hash[:gene_ids]
all_motifs_found = true
motifs_to_include.each {|x|
  if !motifs.keys.include?(x)
    STDERR.puts "#{x} not found in motif definitions."
    all_motifs_found = false
  end
}
exit unless all_motifs_found


status_string = "Building model composed of unbound"
status_string << ", nucleosome" if args_hash[:include_nucleosome]
status_string << ", #{motifs_to_include.join(', ')}" if motifs_to_include.length > 0
status_string << "."
STDERR.puts status_string

# array of transition probabilities from distributor to background, then to each motif
initial_transitions = [1.0]
motifs_to_include.length.times { initial_transitions.push(0.01) }
initial_transitions.push(1.0) if args_hash[:include_nucleosome]

n_emitting_states = 1
motifs_to_include.each {|m| n_emitting_states += 2 * motifs[m].length}

e = read_dinucleotide_freqs("in_vitro_nucleosome_dinucleotide_probs.txt")
if args_hash[:include_nucleosome]
  foo = construct_nucleome_model(e, background_motif, model_alphabet, n_emitting_states - 1, initial_transitions[-1])
  n_emitting_states += foo[0]
end

n_states = n_emitting_states + 1 + motifs_to_include.length
# background + distributor + (2*motif_length + 1 silent state) per motif
# background will always be state 0, then first motif is 1 through its length
# distributor will always be state n_emitting_states


puts "model = {"
puts "  n_states = #{n_states};"
puts "  silent_states_begin = #{n_emitting_states};"
puts "  alphabet_length = #{model_alphabet.length};"
puts "  alphabet = \"#{model_alphabet.join}\";"

# use uniform probability of beginning in any emitting state, for now
puts "  initial_probs = ("
p = 1.0 / n_emitting_states
(0..n_emitting_states-2).each {|i|
  puts "                   (#{i}, #{p}),"
}
puts "                   (#{n_emitting_states-1}, #{p})"
puts "  );"

# output state transition matrix

index = 0
motif_starts = []
puts ""
puts "  transition_matrix = ("

# non-silent states
puts "                       // non-silent states"
puts "                       // background to distributor"
puts "                       (#{index}, #{n_emitting_states}, 1.0),"
index += 1
motifs_to_include.each {|m|
  motif_starts.push(index)
  forward = true
  2.times {
    if (forward)
      puts "                       // traverse along motif for #{m} (forward)"
      forward = false
    else
      puts "                       // traverse along motif for #{m} (reverse)"
    end
    (motifs[m].length - 1).times {
      puts "                       (#{index}, #{index + 1}, 1.0),"
      index += 1
    }
    puts "                       (#{index}, #{n_emitting_states}, 1.0),"
    index += 1
  }
}

puts foo[1] if args_hash[:include_nucleosome]

# silent states
puts "                       // silent states"
puts "                       // distributor to background"
puts "                       (#{n_emitting_states}, 0, #{initial_transitions[0]})#{motifs_to_include.length > 0 ? ',' : ''}"

(1..(motifs_to_include.length)).each {|c|
  m_name = motifs_to_include[c-1]
  puts "                       // distributor to #{m_name}"
  puts "                       (#{n_emitting_states}, #{n_emitting_states + c}, #{initial_transitions[c]}),"
  puts "                       // #{m_name} silent state to forward motif"
  puts "                       (#{n_emitting_states + c}, #{motif_starts[c-1]}, 0.5),"
  puts "                       // #{m_name} silent state to reverse motif"
  # account for no comma on the last entry...
  comma = ","
  comma = "" if c == motifs_to_include.length
  puts "                       (#{n_emitting_states + c}, #{motif_starts[c-1] + motifs[m_name].length}, 0.5)#{comma}"
}

puts "                      );"



# output emission matrix
index = 0
puts ""
puts "  emission_matrix = ("

puts foo[2] if args_hash[:include_nucleosome]

# output background emissions
puts "                     // background"
model_alphabet.each_index {|i|
  puts "                     (0, #{i}, #{background_motif[i]})#{(i < model_alphabet.length - 1 || (motifs_to_include.length > 0)) ? ',' : ''}"
}

index += 1
motifs_to_include.each {|m|
  motif_starts.push(index)
  puts "                     // traverse along motif for #{m} (forward)"
  motifs[m].each_index {|m_pos|
    model_alphabet.each_index {|i|
      motif_alphabet_i = motifs_alphabet_map[model_alphabet[i]]
      puts "                     (#{index}, #{i}, #{motifs[m][m_pos][motif_alphabet_i]}),"
    }
    index += 1
  }

  puts "                     // traverse along motif for #{m} (reverse)"
  (0..(motifs[m].length-1)).to_a.reverse.each {|m_pos|
    model_alphabet.each_index {|i|
      motif_alphabet_i = motifs_alphabet_map[reverse_map[model_alphabet[i]]]
      # account for no comma on the last entry...
      comma = ","
      comma = "" if (m == motifs_to_include[-1] && m_pos == 0 && i == model_alphabet.length - 1)
      puts "                     (#{index}, #{i}, #{motifs[m][m_pos][motif_alphabet_i]})#{comma}"
    }
    index += 1
  }
}

puts "                    );"


puts "};"
