import json


# capture different options 
input_dict = {
    "intermediate_structure" : "ortho", # options are ortho, para, meta
    "number_of_theta_values" : "5", # integers > 0
    "number_of_phi_values" : "5", # integers > 0
    "orbital_basis" : "6-311+G*",
    "omega_value_ev" : 1.8,
    "lambda_magnitude_au" : 0.1,
    "lambda_vector" : [1.0, 0, 0],
    "para_coordinates" :  [
 "C                 -0.80658313    1.22973465    0.03041801",
 "C                  0.56153576    1.23725234    0.01622618",
 "C                  1.22915389    0.01001055    0.01220575",
 "H                 -1.36676923    2.15803094    0.04420367",
 "H                  1.14116413    2.14927050    0.01037697",
 "N                  2.71357475    0.03144573   -0.00289824",
 "O                  3.28013247   -1.09741954   -0.00254733",
 "O                  3.24714953    1.17621948   -0.01252002",
 "C                 -0.77042978   -1.26805414    0.04039660",
 "H                 -1.30353926   -2.21202933    0.06122375",
 "C                  0.59726287   -1.23605918    0.02634378",
 "H                  1.20308359   -2.13089607    0.02793117",
 "C                 -1.56287141   -0.03049318    0.01040538",
 "H                 -2.41148563   -0.03994459    0.70143946",
 "Br                -2.40993182   -0.04931830   -1.82359612"],
 "meta_coordinates" : [
"C                  0.02949981    1.33972592    0.06817723",
 "C                  1.43483278    1.28667967    0.00635313",
 "C                  2.11179024    0.05106117   -0.00544138",
 "C                  1.44506636   -1.13720058    0.03116583",
 "C                 -0.68793171    0.16822220    0.10995314",
 "H                 -0.47126997    2.29839666    0.07811355",
 "H                  2.02732783    2.19651728   -0.03220624",
 "H                  1.98966526   -2.07643217    0.02318494",
 "H                 -1.77163480    0.18040547    0.15819632",
 "N                  3.58635895    0.05097292   -0.06745286",
 "O                  4.14711759   -1.05966097   -0.08807849",
 "O                  4.14497859    1.16390951   -0.09010823",
 "C                 -0.02361177   -1.14582791    0.08353483",
 "H                 -0.43674996   -1.87247364    0.78889576",
 "Br                -0.53591638   -1.86972195   -1.74078671"
 ],
 "ortho_coordinates" : [
 "C                  0.51932475    1.23303451   -0.03194925",
 "C                  1.94454413    1.26916358   -0.03672882",
 "C                  2.62037793    0.09283428   -0.02499003",
 "C                 -0.19603352    0.03013062    0.00102732",
 "H                 -0.02069420    2.17423764   -0.04336646",
 "H                  2.48281698    2.20891057   -0.03611879",
 "H                 -1.27770137    0.03990295    0.01166953",
 "N                  4.09213475    0.09594076    0.03662979",
 "O                  4.63930696   -1.02169275    0.14459220",
 "O                  4.66489883    1.19839699   -0.02327545",
 "C                  0.49428518   -1.16712649    0.02099746",
 "H                 -0.03251071   -2.11492669    0.05447935",
 "C                  1.96291176   -1.21653219   -0.02111314",
 "H                  2.44359113   -1.96306433    0.61513886",
 "Br                 2.17304025   -1.94912156   -1.90618750"
 ],
}




def write_json(filename, data):
    with open(filename, 'w') as f:
        json.dump(data, f, indent=4)


write_json('input_options.json', input_dict)