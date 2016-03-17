

from os.path import expanduser
from visualization import read_occupancy_profile, plot_occupancy_profile

if __name__ == '__main__':

    op = read_occupancy_profile('~/COMPETE_examples/YIR022W_IX_397182_399867.tsv')
    chromo = 9
    coordinate_start = 397182
    
    macisaac_file = '~/COMPETE_examples/macisaac_p005_c2.csv'
    orf_file = '~/COMPETE_examples/saccharomyces_cerevisiae.20080621.gff'

    plot_occupancy_profile(op, chromo = chromo, 
        coordinate_start = coordinate_start, 
        # orf_annotation and macisaac_annotation can be set to None if not needed
        orf_annotation = orf_file, macisaac_annotation = macisaac_file, 
        file_name = expanduser('~/COMPETE_examples/YIR022W_IX_397182_399867.tsv.png'))
