#!/usr/bin/env ruby

require 'yaml'
require 'pp'

config = YAML.load(ARGF)
expr_means = YAML.load_file("gene_expr_means.txt")
expr_means.delete_if {|k, v| v.class != Float}
baseline = expr_means.values.min
expr_means.each {|k, v| expr_means[k] = v - baseline}

expr_means.keys.sort.each {|k|
  k = k.downcase
  config['motifs'] << k
  config['conc'][k] = 0.0
}
config['motifs'].uniq!

config['conc'].each_key{|k| config['conc'][k] = (2**expr_means[k.upcase]) if expr_means[k.upcase]}

YAML.dump(config, STDOUT)
