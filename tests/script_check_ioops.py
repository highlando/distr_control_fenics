from dolfin_navier_scipy.problem_setups import cyl_fems
import optical_test_outpop as too
import optical_test_inpop as tio

N = 3
NU = 5
NY = 5
femp = cyl_fems(N)

too.check_output_opa(NYx=1, NYy=2, femp=femp)
tio.check_input_opa(NU, femp)
