###### Ignore this
import pytest
import hashlib
_ = 123456789  # just a wrong number, please ignore
###### Stop ignoring

# Some stuff you can/should use ...
# feel free to import additional things from those packages already imported
# or the Python Standard Library (https://docs.python.org/3/library/)
# (if it helps) but do not import other 3rd party packages.

import numpy as np
from cobra.io import read_sbml_model
from cobra import Reaction, Metabolite
from cobra.util import create_stoichiometric_matrix

# Read model (central metabolism model of Escherichia coli); don't change
model = read_sbml_model("e_coli_core.xml")
medium = model.medium
del(medium["EX_co2_e"])
model.medium = medium

# General hints:
# 1. Use the E. coli Core model (`model`) in its default configuration as configured above if not stated otherwise.
# 2. Remember to undo modifications to the model before continuing with the next task
#    (either make a copy of the model for each task or use the `with model: ...` construct as shown in the exercise).


# 1. Based on the model, what is the theoretical maximum yield of acetate in units of
#    mmol-acetate/mmol-glucose?
# Hints:
# * Needs to be a positive floating point number

# Put your intermediate solution steps here if you have any ...

# Replace _ with you're final calculation step or a variable that contains the final solution.
# maximum_acetate_yield_mol needs to resolve to a positive floating point number
maximum_acetate_yield_mol = _


# 2. Based on the model, what is the theoretical maximum yield of acetate in units of
#    cmol-acetate/cmol-glucose?
# Hints:
# * You can look up the elemental composition of a metabolite usinge `.elements` (which is a dict)

# Put your intermediate solution steps here if you have any ...


# Replace _ with you're final calculation step or a variable that contains the final solution.
# maximum_acetate_yield_cmol needs to resolve to a positive floating point number
maximum_acetate_yield_cmol = _


# 3. Based on the model's stoichiometry alone, how many reaction fluxes need to be measured
#    in order to make the system determined and solvable?
# Hints:
# * You can use `create_stoichiometry matrix` (imported above) to extract a stoichiometric matrix
#   from the model.
# * You can use numpy (imported as `np` above) in case you need some matrix related functionality.

# Put your intermediate solution steps here if you have any ...


# Replace _ with you're final calculation step or a variable that contains the final solution.
# degrees_of_freedom needs to be an integer number.

degrees_of_freedom = _


# 4. How much is the (optimal) growth rate reduced if fumarase (FUM in the model) is
#    overexpressed to have a 2-fold higher flux in comparison to its flux at maximum growth rate?

# Hints:
# * 2-fold means the reference flux of FUM (at optimal growth) multiplied by two
# * Growth reduction should be calculated as optimal_growth - growth_fum_overexpression

# Put your intermediate solution steps here if you have any ...


# Replace _ with you're final calculation step or a variable that contains the final solution.
# growth_reduction needs to be a floating point number that corresponds to
# optimal_growth - growth_fum_overexpression

growth_reduction = _


# 5. What genes are essential under acetate conditions but not glucose conditions?
# Hints:
# * The biomass objective already set in the model should be used to determine the maximum growth rate
# * For the glucose condition use the models as is; for the acetate condition remove
#   glucose ('EX_glc__D_e') from the medium and set a 10 for 'EX_ac_e' (corresponds to a maximum uptake
#   rate of 10 mmol gDW^-1 h^-1)
# * Use `model.slim_optimize(error_value=0.)` to determine the mutant growth rate
# * You can knock out a gene using its `.knock_out()` method
# * Assume a gene is essential if the growth rate drops < 0.05 h^-1 upon being knocked out
# * Remember to use `with model: ...` to revert changes like knockouts
# * You can use Python sets (https://docs.python.org/3.8/library/stdtypes.html#set-types-set-frozenset);
#   in particular the `difference` method should be useful, e.g.,
#   set(['a', 'b', 'c']).difference(set['c', 'd']) ==> set(['a', 'b'])

# Put your intermediate solution steps here if you have any ...

    
# Replace _ with you're final calculation step or a variable that contains the final solution.
# essential_only_in_acetate needs to resolve to a set of gene IDs of type str ({'b2286', 'b2287', ...})
essential_only_in_acetate = _


#### Tests are happening in the end now ...
###### Don't touch

def test_maximum_acetate_yield_mol():
    assert maximum_acetate_yield_mol == pytest.approx(2.)
    
def test_maximum_acetate_yield_cmol():
    assert maximum_acetate_yield_cmol == pytest.approx(0.6666666666666666)

def test_degrees_of_freedom():
    assert degrees_of_freedom == 28

def test_growth_reduction():
    assert growth_reduction == pytest.approx(0.2244331576552907)

def test_essential_only_in_acetate():
    assert essential_only_in_acetate == {'b2286', 'b2287', 'b3737', 'b2282', 'b2281', 'b2283', 'b3731', 'b2279', 'b3919', 'b3735', 'b3734', 'b3736', 'b0722', 's0001', 'b2277', 'b2288', 'b3732', 'b3738', 'b2284', 'b0721', 'b2278', 'b4025', 'b0723', 'b0724', 'b2285', 'b2280', 'b4015', 'b3733', 'b2276'}
###### this
