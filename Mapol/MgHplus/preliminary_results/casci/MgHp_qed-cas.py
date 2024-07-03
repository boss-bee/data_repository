import numpy as np
import sys
sys.path.append("/home/nvu12/software/qed_ci_main/qed_ci_directci2/qed-ci/src")
#np.set_printoptions(precision=5, linewidth=200, suppress=True)
import psi4
import json
from helper_PFCI import PFHamiltonianGenerator
from helper_PFCI import Determinant
from helper_cqed_rhf import cqed_rhf
np.set_printoptions(threshold=sys.maxsize)

mol_tmpl = """
    1 1
    Mg 0.000 0.0000 0.000
    H  0.000 0.0000 **R**
    symmetry c1
"""
mol_str = """
    1 1
    Mg 0.000 0.0000 0.000
    H  0.000 0.0000 **R**
    symmetry c1
"""
options_dict = {'basis': 'cc-pVTZ',
                  'scf_type': 'pk',
                  'e_convergence': 1e-10,
                  'd_convergence': 1e-10
                  }

r_array = np.linspace(1.25, 4.0, 50)

cavity_dict = {
    'omega_value' : 0.,
    'lambda_vector' : np.array([0., 0., 0.]),
    #'omega_value' : 0.16537186620313543,
    #'lambda_vector' : np.array([0, 0, 0.05]),
    'ci_level' : 'cas',
    'davidson_roots' : 10,
    'davidson_maxdim':6,
    'nact_orbs' : 10,
    'nact_els' : 12,
    'photon_number_basis' : False,
    'canonical_mos' : False,
    'coherent_state_basis' : True,
    'full_diagonalization' : False,
    # specify the number of photons - 2 means |0> , |1>, |2> will be in the basis
    'number_of_photons' : 0
}

# set up base dictionary - some of this will be updated with each calculation
dictionary = {
  "repo" : {
      
      "repository_url" : "https://github.com/mapol-chem/qed-ci",
      "branch" : "direct_ci",
      "commit" : "7f25d1c5fd458a06002e79d610aa17399f213f71"
  },
    
  "molecule" : {
    "molecule_name": "MgH",
    "geometry": [],
    "symbols": [
      "Mg",
      "H"
    ],
    "bond_length": [],
  },
    "driver": "energy",
    "model" : {
        "method" : "qed-cas",
        "orbital_basis"  : "cc-pVTZ",
        "photon_basis" : "coherent_state_basis",
        "number_of_photon_states" : cavity_dict["number_of_photons"],
        "nact_orbs" : cavity_dict["nact_orbs"],
        "nact_els" : cavity_dict["nact_els"],
        "lambda" : [
            cavity_dict["lambda_vector"][0],
            cavity_dict["lambda_vector"][1],
            cavity_dict["lambda_vector"][2]
        ],
        "omega" : cavity_dict["omega_value"]  
    },
    
    "return_result" : [
        
    ],
    
}

# function to generate file names based on system details
def generate_file_name(dic):
    file_name = dic["molecule"]["molecule_name"] + "_"
    file_name += str(dictionary["model"]["method"]) + "_"
    file_name += str(dictionary["model"]["orbital_basis"]) + "_"
    file_name += str(dictionary["model"]["photon_basis"]) + "_"
    file_name += str(dictionary["model"]["number_of_photon_states"]) + "_"
    file_name += str(dictionary["model"]["nact_els"]) + "_"
    file_name += str(dictionary["model"]["nact_orbs"]) + "_"
    file_name += "lambda_" + str(dictionary["model"]["lambda"][0]) 
    file_name += "_" + str(dictionary["model"]["lambda"][1]) 
    file_name += "_" + str(dictionary["model"]["lambda"][2]) + "_"
    file_name += "omega_" + str(dictionary["model"]["omega"]) + ".json"
    return file_name

psi4.set_options(options_dict)
psi4.core.set_output_file('mgh_1.400.out', True)

#LiH = PFHamiltonianGenerator(mol_str, options_dict, cavity_dict)
#print(np.shape(LiH.H_PF))
#for i in range(10):
# print('state', i, LiH.CIeigs[i])
file_name = generate_file_name(dictionary) 
for r in r_array:
    mol_str = mol_tmpl.replace("**R**", str(r))
    mol = psi4.geometry(mol_str)
    LiH = PFHamiltonianGenerator(mol_str, options_dict, cavity_dict)
    dictionary["return_result"].append(list(LiH.CIeigs))
    #geometry = np.asarray(mol.geometry())
    #geometry_list = geometry.tolist()
    #dictionary["molecule"]["geometry"].append(geometry_list)
    dictionary["molecule"]["bond_length"].append(r)

print(dictionary)

json_object = json.dumps(dictionary, indent=4)

with open(file_name, "w") as outfile:
    outfile.write(json_object)
