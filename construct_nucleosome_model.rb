#!/usr/bin/env ruby

class Array
  def sum
    inject(0) { |s, v| s += v }
  end
  def /(d)
    collect {|i| i.to_f/d}
  end
end

def construct_nucleome_model(emissions, background_motif, model_alphabet, start_index, p_trans_in)
  index = start_index + 1
  trans = ""
  emit = ""
  # emissions is an array of hashes with keys being dinucleotides and values being their probability of emission

  n_bg_padding = (147-emissions.length)/2 + 5
#  n_bg_padding = 8

  trans << "                       // begin nucleosome model\n"
  emit << "                     // begin nucleosome model\n"

  # output n_bg_padding background states
  trans << "                       // #{n_bg_padding-1} normal background states\n"
  emit << "                     // #{n_bg_padding-1} normal background states\n"
  model_alphabet.each_index {|i|
    emit << "                     (#{index}, #{i}, #{background_motif[i]}),\n"
  }
  index += 1
  (n_bg_padding-2).times {
    trans << "                       (#{index-1}, #{index}, 1.0),\n"
    model_alphabet.each_index {|i|
      emit << "                     (#{index}, #{i}, #{background_motif[i]}),\n"
    }
    index += 1
  }

  trans << "                       // 1 branched background state\n"
  emit << "                     // 1 branched background state\n"
  model_alphabet.each_index {|i|
    trans << "                       (#{index-1}, #{index+i}, #{background_motif[i]}),\n"
    emit << "                     (#{index+i}, #{i}, 1.0),\n"
  }
#  index += model_alphabet.length

  pos = 0 # nucleosome position
  model_alphabet.each_index {|i|
    model_alphabet.each_index {|j|
      dest_index = index + model_alphabet.length * (1 + i) + j
      trans << "                       (#{index+i}, #{dest_index}, #{emissions[pos][i][j]}),\n"
      emit << "                     (#{dest_index}, #{j}, 1.0),\n"
    }
  }
  index += model_alphabet.length * (1 + model_alphabet.length)
  pos += 1


  # output (emissions.length - 1) * 16 nucleosome states
  (emissions.length - 1).times {
    model_alphabet.each_index {|i|
      model_alphabet.each_index {|j|
        src_index = index - (model_alphabet.length ** 2) + i * model_alphabet.length + j
        emitting_state = index + i * model_alphabet.length + j
        emit << "                     (#{emitting_state}, #{j}, 1.0),\n"
        model_alphabet.each_index {|k|
          dest_index = index + model_alphabet.length * j + k
          trans << "                       (#{src_index}, #{dest_index}, #{emissions[pos][j][k]}),\n"
        }
      }
    }
    index += model_alphabet.length ** 2
    pos += 1
  }

  # make all 16 of the last nucleosome states go to next background state
  old_index = index - model_alphabet.length ** 2
  (model_alphabet.length ** 2).times {
    trans << "                       (#{old_index}, #{index}, 1.0),\n"
    old_index += 1
  }

  # output n_bg_padding background states
  trans << "                       // #{n_bg_padding} normal background states\n"
  emit << "                     // #{n_bg_padding} normal background states\n"
  model_alphabet.each_index {|i|
    emit << "                     (#{index}, #{i}, #{background_motif[i]}),\n"
  }
  index += 1
  (n_bg_padding-1).times {
    trans << "                       (#{index-1}, #{index}, 1.0),\n"
    model_alphabet.each_index {|i|
      emit << "                     (#{index}, #{i}, #{background_motif[i]}),\n"
    }
    index += 1
  }

  n_nuc_states = index - start_index - 1
  trans << "                       // distributor to nucleosome\n"
  trans << "                       (#{start_index+n_nuc_states+1}, #{start_index+1}, #{p_trans_in}),\n"
  trans << "                       // nucleosome to distributor\n"
  trans << "                       (#{index-1}, #{start_index+n_nuc_states+1}, 1.0),\n"
  [n_nuc_states, trans, emit]
end


def read_dinucleotide_freqs(filename)
  d = Array.new
  a = Array.new
  File.new(filename).readlines.each{|l| d.push(l.chomp.split(/\t/).collect{|i| i.to_f})}

  a = d.collect {|l| 
    s1 = l[0..3].sum
    s2 = l[4..7].sum
    s3 = l[8..11].sum
    s4 = l[12..15].sum

    [l[0..3]/s1, l[4..7]/s2, l[8..11]/s3, l[12..15]/s4]
  }

  a
end
