import numpy as np
import psi4
import json
from helper_PFCI import PFHamiltonianGenerator
import sys
np.set_printoptions(threshold=sys.maxsize)


mol_str = """
0 1
C            0.000004551444    -0.711108766560     0.000010757138
C           -0.000004765464     0.711109074728     0.000001509641
C            1.231111212344    -1.392932324772    -0.000008872321
C            1.231092818353     1.392949010021    -0.000005886632
C           -1.231092617348    -1.392949084219     0.000002684292
C           -1.231111293835     1.392932340511     0.000007560162
C            2.416252154882    -0.702767074233    -0.000017781801
C            2.416242565430     0.702799552187    -0.000000723902
C           -2.416242703804    -0.702799153654     0.000010510747
C           -2.416251802225     0.702766484279     0.000002457080
H            1.229116150588    -2.471858679894    -0.000015215757
H            1.229083140542     2.471874840369    -0.000007074480
H           -1.229084430358    -2.471875468868    -0.000008700875
H           -1.229118355038     2.471858460776     0.000007352885
H            3.350158261997    -1.238508894745    -0.000013806252
H            3.350141246250     1.238554060765     0.000000769827
H           -3.350140729721    -1.238555163085     0.000018452996
H           -3.350156710411     1.238510150658    -0.000008144872
symmetry c1
"""


options_dict = {'basis': 'cc-pVDZ',
                  'scf_type': 'pk',
                  'e_convergence': 1e-11,
                  'd_convergence': 1e-11
                  }

l_array = np.array([0., 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01])

cavity_dict = {
    'omega_value' : 0.21767223258098056,
    'lambda_vector' : np.array([0., 0.0, 0.0]),
    'ci_level' : 'cas',
    'davidson_roots' : 15,
    'davidson_maxdim':6,
    'photon_number_basis' : False,
    'canonical_mos' : False,
    'coherent_state_basis' : True,
    'nact_orbs' : 10,
    'nact_els' : 10,
    'full_diagonalization' : False,
    'number_of_photons' : 5
}

# set up base dictionary - some of this will be updated with each calculation
dictionary = {
  "repo" : {
      
      "repository_url" : "https://github.com/mapol-chem/qed-ci",
      "branch" : "direct_ci",
      "commit" : "7f25d1c5fd458a06002e79d610aa17399f213f71"
  },
    
  "molecule" : {
    "molecule_name": "C10H8",
    "geometry": [mol_str],
    "symbols": [
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H"
    ],
    "coordinates": "\n0 1\nC            0.000004551444    -0.711108766560     0.000010757138\nC           -0.000004765464     0.711109074728     0.000001509641\nC            1.231111212344    -1.392932324772    -0.000008872321\nC            1.231092818353     1.392949010021    -0.000005886632\nC           -1.231092617348    -1.392949084219     0.000002684292\nC           -1.231111293835     1.392932340511     0.000007560162\nC            2.416252154882    -0.702767074233    -0.000017781801\nC            2.416242565430     0.702799552187    -0.000000723902\nC           -2.416242703804    -0.702799153654     0.000010510747\nC           -2.416251802225     0.702766484279     0.000002457080\nH            1.229116150588    -2.471858679894    -0.000015215757\nH            1.229083140542     2.471874840369    -0.000007074480\nH           -1.229084430358    -2.471875468868    -0.000008700875\nH           -1.229118355038     2.471858460776     0.000007352885\nH            3.350158261997    -1.238508894745    -0.000013806252\nH            3.350141246250     1.238554060765     0.000000769827\nH           -3.350140729721    -1.238555163085     0.000018452996\nH           -3.350156710411     1.238510150658    -0.000008144872\nsymmetry c1"
    },
    "driver": "energy",
    "model" : {
        "method" : "qed-cas",
        "orbital_basis"  : options_dict["basis"],
        "photon_basis" : "coherent_state_basis",
        "number_of_photon_states" : cavity_dict["number_of_photons"],
        "nact_orbs" : 10,
        "nact_els" : 10,
        "lambda" : [
        ],
        "omega" : cavity_dict["omega_value"]  
    },
    
    "return_result" : [
        
    ],
    
}

# function to generate file names based on system details
def generate_file_name(dictionary):
    file_name = dictionary["molecule"]["molecule_name"] + "_"
    file_name += str(dictionary["model"]["method"]) + "_"
    file_name += str(dictionary["model"]["orbital_basis"]) + "_"
    file_name += str(dictionary["model"]["photon_basis"]) + "_"
    file_name += str(dictionary["model"]["number_of_photon_states"]) + "_"
    file_name += str(dictionary["model"]["nact_els"]) + "_"
    file_name += str(dictionary["model"]["nact_orbs"]) + "_"
    #file_name += "lambda_" + str(dictionary["model"]["lambda"][0]) 
    #file_name += "_" + str(dictionary["model"]["lambda"][1]) 
    #file_name += "_" + str(dictionary["model"]["lambda"][2]) + "_"
    file_name += "omega_" + str(dictionary["model"]["omega"]) + ".json"
    return file_name

psi4.set_options(options_dict)
psi4.core.set_output_file('h2o_1.400.out', True)



file_name = generate_file_name(dictionary) 
mol = psi4.geometry(mol_str)
#LiH = PFHamiltonianGenerator(mol_str, options_dict, cavity_dict)
#dictionary["return_result"].append(list(LiH.CIeigs))

for l in l_array:
    l_vector = np.array([0., l, 0])
    cavity_dict['lambda_vector'] = np.copy(l_vector)
    lam_list = [0, l, 0]
    dictionary["model"]["lambda"].append(lam_list)
    mol = psi4.geometry(mol_str)
    LiH = PFHamiltonianGenerator(mol_str, options_dict, cavity_dict)
    dictionary["return_result"].append(list(LiH.CIeigs))

print(dictionary)

json_object = json.dumps(dictionary, indent=4)

with open(file_name, "w") as outfile:
    outfile.write(json_object)
