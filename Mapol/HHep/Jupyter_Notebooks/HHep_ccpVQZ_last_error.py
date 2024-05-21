import numpy as np
import sys
import psi4
from helper_PFCI import PFHamiltonianGenerator
np.set_printoptions(threshold=sys.maxsize)
psi4.core.set_output_file('output.dat', False)
import time
import json

qed_fci_qz_large_lambda = np.array([
 -2.9773960939e+00,
 -2.1846350116e+00,
 -2.0247781044e+00,
 -2.0023579319e+00,
 -1.7700528337e+00,
 -1.7270874847e+00,
 -1.7270874847e+00,
 -1.7220652946e+00,
 -1.6602926834e+00,
 -1.6602926834e+00]
)





# these file names should still be good
qz_en_file = "HHep_fci_cc-pVQZ_Energies.npy"
qz_mu_file = "HHep_fci_cc-pVQZ_Dipoles.npy"


qz_en = np.load(qz_en_file)
qz_mu = np.load(qz_mu_file)

# setup basic arguments to create an instance of the PFHamiltonianGenerator class
mol_str = """
Li
H 1 1.4
symmetry c1
"""

options_dict = {
    "basis": "sto-3g",
    "scf_type": "pk",
    "e_convergence": 1e-10,
    "d_convergence": 1e-10,
}


cavity_free_dict = {
    'omega_value' : 0.0,
    'lambda_vector' : np.array([0, 0, 0.0]),
    'ci_level' : 'fci',   
    'full_diagonalization' : True,
    'number_of_photons' : 0, 
}

# create the instance of our PFHamiltonianGenerator class
qz_inst = PFHamiltonianGenerator(mol_str, options_dict, cavity_free_dict)

# number of adiabatic states for each basis set
N_el_qz = len(qz_en)

# define 95% and 100% of electronic states
N_el_95_pct = 3420
N_el_100_pct = 3600

# number of photonic states
N_ph = 3

# omega
omega_qz = 0.9637811053

# lambda vector
lambda_vector = np.array([0., 0., 0.02])


# build PCQED Hamiltonian
qz_inst.fast_build_pcqed_pf_hamiltonian(N_el_100_pct, N_ph, omega_qz, lambda_vector, qz_en, qz_mu, neglect_DSE=False)

ei = np.copy(qz_inst.PCQED_pf_eigs[:6])
print(F' {N_el_100_pct}, {ei[0]:12.10e},{ei[1]:12.10e},{ei[2]:12.10e},{ei[3]:12.10e},{ei[4]:12.10e},{ei[5]:12.10e}')
