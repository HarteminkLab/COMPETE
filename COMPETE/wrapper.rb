#!/usr/bin/env ruby

require 'yaml'
#require 'pp'

config = YAML.load(ARGF)
#pp config

if config['conc_scale']
  config['conc'].each{|k, v|
    config['conc'][k] = v * config['conc_scale']
  }
end

construct_model = "./construct_model_from_motifs.rb #{config['motifs'].join(' ')}"
construct_model << " -n" if config['nuc']
construct_model << " > #{config['output_basename']}_model.cfg"
STDERR.puts construct_model
system(construct_model)

f = File.new("#{config['output_basename']}_seq_filenames.txt", "w")
f.puts "chr/#{config['chr'].upcase}.txt #{config['range'][0]} #{config['range'][1]}"
f.close
STDERR.puts "#{config['output_basename']}_seq_filenames.txt: chr/#{config['chr'].upcase}.txt #{config['range'][0]} #{config['range'][1]}"


run_model = "./compete -t #{config['inverse_temp']} -u #{config['conc']['unbound']}"
output_filename = "#{config['output_basename']}_chr#{config['chr']}_#{config['range'][0]}..#{config['range'][1]}_unbound_#{config['conc']['unbound']}"

if config['nuc']
  run_model << " -n #{config['conc']['nuc']}"
  output_filename << "_nuc_#{config['conc']['nuc']}"
end

if config['motifs'].length > 0
  motif_conc_param = config['motifs'].collect{|x| config['conc'][x]}.join(',')
  run_model << " -m #{motif_conc_param}"
  run_model << " -N #{config['motifs'].join(',')}"

  motif_conc_substr = config['motifs'].collect{|x| "_#{x}_#{config['conc'][x]}"}.join('')
  output_filename << motif_conc_substr unless motif_conc_substr.length >= 64
end

if config['fix_positions'] && config['fix_positions'].length > 0
  fix = Array.new
  motif_states = Hash.new
  motif_states["unbound"] = 0

  d = IO.readlines("#{config['output_basename']}_model.cfg").join
  if config['nuc']
    md = d.match(/distributor to nucleosome.*?\(\d+, (\d+)/m)
    motif_states["nuc"] = md[1].to_i
  end
  if config['motifs'].length > 0
    config['motifs'].each{|m|
      md = d.match(/#{m.upcase} silent state to forward motif.*?\(\d+, (\d+)/m)
      motif_states[m] = md[1].to_i
    }
  end

  config['fix_positions'].each{|k, v|
    fix.push("#{k[0]-config['range'][0]+1}-#{k[1]-config['range'][0]+1}:#{motif_states[v]}")
  }

  run_model << " -f #{fix.join(',')}"
end

if config['output_start_probs_only']
  run_model << " -s"
end

run_model << " #{config['output_basename']}_model.cfg #{config['output_basename']}_seq_filenames.txt"
output_filename << ".txt"
output_filename = config['output_filename'] if config['output_filename']
run_model << " > #{output_filename}"
STDERR.puts
STDERR.puts run_model
system(run_model)
