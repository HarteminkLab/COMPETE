
'''
define color related code
'''

colorset48 = ['#8c2a2a', '#a69863', '#26803e', '#0066ff', '#660052', 
            '#f22000', '#bfb639', '#49f2a3', '#263e80', '#ff4dc4', '#401100',
            '#eeff00', '#4dfff3', '#0022ff', '#bf73a6', '#bf5d39', '#535900',
            '#2e9992', '#e6e9ff', '#b20047', '#ff944d', '#d4d9c3', '#0f3133',
            '#676973', '#bfacb4', '#ffc299', '#ddff99', '#e6fdff', '#000059',
            '#ff0044', '#733d00', '#66bf00', '#91dff2', '#6836b3', '#591b2b',
            '#403d39', '#748073', '#0099e6', '#a582d9', '#e68a96', '#d99100',
            '#00ff00', '#005580', '#3a2e4d', '#594d36', '#134019', '#001b33', '#d600e6']

predefined_dbfs = ['Abf1', 'Ace2', 'Adr1', 'Aft2', 'Cbf1', 'Cin5', 'Cst6', 'Dal80', 
                    'Dal82', 'Fhl1', 'Fkh1', 'Gat1', 'Gcn4', 'Gcr1', 'Gln3', 'Hap1', 
                    'Hcm1', 'Mbp1', 'Mcm1', 'Met31', 'Mot3', 'Msn4', 'Nrg1', 'Phd1', 
                    'Pho2', 'Rap1', 'Reb1', 'Rox1', 'Sfp1', 'Skn7', 'Smp1', 'Sok2', 
                    'Ste12', 'Sut1', 'Swi4', 'Swi5', 'Tec1', 'Ume6', 'Xbp1', 'Yap6', 
                    'Yap7', 'Yox1', 'Gal4', 'Pho4', 'Arr1', 'Mac1', 'Ndd1', 'Yap5']

nucleosome_color = '0.7'

default_dbf_color_map = dict(zip(predefined_dbfs, colorset48))
default_dbf_color_map['nucleosome'] = nucleosome_color