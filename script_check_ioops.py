from dolfin_navier_scipy.problem_setups import cyl_fems
import distr_control_fenics.optical_test_outpop as too
import distr_control_fenics.optical_test_inpop as tio

N = 3
NU = 5
NY = 5
femp = cyl_fems(N)

# too.check_output_opa(NY, femp)
tio.check_input_opa(NU, femp)
