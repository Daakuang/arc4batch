"""Retain the validated physical model/EKF; align the optimization with temperature quality."""
import casadi as ca
from parameter_model import ParameterModel


class QualityModel(ParameterModel):
    def __init__(self,temperature_weight=3e6/360):
        super().__init__()
        z=ca.SX.sym('quality_state',12);u=ca.SX.sym('quality_input',3);theta=ca.SX.sym('quality_theta',2)
        dx,oldcost=self.f(z,u,theta)
        # Raise only temperature priority; preserve production and pressure terms.
        # The common 1/1000 scale also applies to move regularization, but not
        # to the strong explicit feasibility-recovery penalty.
        self.f=ca.Function('quality_priority_rhs',[z,u,theta],[dx,oldcost/1000.])
        self.hold_penalty=ca.Function('quality_hold_penalty',[z],
            [(temperature_weight-3e2/360)*(z[5]*self.scale[5]-351.)**2/1000.])
        self.temperature_weight=temperature_weight
