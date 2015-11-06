def parse_macisaac_motif_file(filename)
  h = Hash.new
  d = IO.readlines(filename).join

  calculated_bg = [0.3103490, 0.1907995, 0.3093810, 0.1898473]

  # specially crafted magical regex to parse the specific format of this file
  d.scan(/Log-odds.*?\n.*?\n(.*?)\s*\n(.*?)\s*\n(.*?)\s*\n(.*?)\s*\n.*?Source:\s+(.*?)\s/m).each {|i|
    rows = i[0..3].collect{|j| foo = j.split(/\s+/); foo.shift; foo}
    gene = i[-1].upcase
    h[gene] = []
    rows[0].length.times { h[gene].push([]) }
    rows.each_index {|j|
      rows[j].each_index {|k|
        h[gene][k].push(2**(rows[j][k].to_f) * calculated_bg[j])
      }
    }

    # normalize to sum to 1.0
    sums = h[gene].collect{|x| x.inject( nil ) { |sum,y| sum ? sum+y : y }}
    h[gene] = h[gene].collect{|x|
      s = x.inject( nil ) { |sum,y| sum ? sum+y : y }
      x.collect{|y| y / s}
    }
  }

  h
end

