def parse_pbm_motif_directory(directory, discard_header=true)
  h = Hash.new

  Dir.glob("#{directory}/*.pwm").each{|filename|
    File.basename(filename) =~ /(.+)\.pwm/
    gene = $1.upcase
    rows = Array.new
    d = IO.readlines(filename)
    d.shift if discard_header

    d.each{|r|
      foo = r.split(/\s+/)
      foo.shift
      rows.push(foo.collect{|x| x.to_f})
    }

    h[gene] = []
    rows[0].length.times { h[gene].push([]) }
    rows.each_index {|j|
      rows[j].each_index {|k|
	h[gene][k].push(rows[j][k])
      }
    }

    # normalize to sum to 1.0; largely unnecessary for pbm data but oh well
    sums = h[gene].collect{|x| x.inject( nil ) { |sum,y| sum ? sum+y : y }}
    h[gene] = h[gene].collect{|x|
      s = x.inject( nil ) { |sum,y| sum ? sum+y : y }
      x.collect{|y| y / s}
    }
  }

  h
end

