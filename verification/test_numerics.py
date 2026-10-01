"""Checks for batch-completion metrics and reduced-reactor identities."""
import sys,unittest
from pathlib import Path
import numpy as np
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src/arc4batch'))
from analyze_results import metric

class NumericalContracts(unittest.TestCase):
    def fixture(self):
        t=np.arange(3303);x=np.zeros((len(t),19));u=np.zeros((len(t)-1,3))
        x[:,0]=1.;x[:,1]=10.;x[:,2]=200.;x[:,11]=351.;x[0,11]=340
        x[:,13]=1.5e6;x[:,15]=3249.
        x[3300:,11]=[351.,351.8,351.2]
        x[3300:,13]=[1.501e6,1.502e6,1e6];x[-1,15]=3250.
        return t,x,u

    def test_error_windows_exclude_startup_and_pressure_finish(self):
        t,x,u=self.fixture();r=metric(t,x,u,scenario='N',controller='test',seed=0)
        self.assertAlmostEqual(r['T_MAE_K'],1/3)
        self.assertAlmostEqual(r['T_absmax_K'],.8)
        self.assertAlmostEqual(r['P_MAE_kPa'],1.5)
        self.assertAlmostEqual(r['P_absmax_kPa'],2)
        self.assertAlmostEqual(r['T_exceedance_percent'],100/3302)
        self.assertTrue(r['completed'])
        self.assertEqual(r['T_exceedance_seconds'],1)

    def test_trip_never_counts_as_completed_batch(self):
        t,x,u=self.fixture();x[-1,11]=373.16
        r=metric(t,x,u,scenario='F',controller='test',seed=0)
        self.assertFalse(r['completed']);self.assertEqual(r['termination'],'temperature_trip')

    def test_cold_recipe_departure_is_distinct_from_thermal_exceedance(self):
        t,x,u=self.fixture();x[3300:,11]=[350.,351.,351.2]
        r=metric(t,x,u,scenario='F',controller='test',seed=0)
        self.assertTrue(r['completed'])
        self.assertEqual(r['T_exceedance_seconds'],0)
        self.assertFalse(r['conforming_completed'])
        self.assertAlmostEqual(r['T_recipe_band_exceedance_percent'],100/3)
        self.assertAlmostEqual(r['T_below_recipe_percent'],100/3)

    def test_completion_requires_recipe_and_operating_conformity_for_time_comparison(self):
        t,x,u=self.fixture();x[3300:,11]=[351.,351.1,351.2]
        r=metric(t,x,u,scenario='N',controller='test',seed=0)
        self.assertTrue(r['conforming_completed'])
        x[100,13]=1.601e6
        r=metric(t,x,u,scenario='N',controller='test',seed=0)
        self.assertTrue(r['completed']);self.assertFalse(r['conforming_completed'])

    def test_missing_charge_prevents_completion(self):
        t,x,u=self.fixture();x[-1,15]=3249.998
        self.assertFalse(metric(t,x,u,scenario='N',controller='test',seed=0)['completed'])

    def test_exact_first_order_coordinate(self):
        NA,NB,V,u,k,C,cf=sp.symbols('NA NB V u k C cf',positive=True)
        reaction=k*NA*NB/V
        derivative=sp.diff(-C*NB/V,NB)*(cf*u-reaction)+sp.diff(-C*NB/V,V)*u
        expected=C*k*NA*NB/V**2-C*(cf-NB/V)*u/V
        self.assertEqual(sp.simplify(derivative-expected),0)

    def test_joint_error_allowance_retains_both_boundary_signs(self):
        C=.09;k=9.72e4*np.exp(-82500/(8.314*350));cf=2.5/.00023
        kp=5e-7;rp=2.777777777777778e-8;rm=.3*rp
        low=(397-350-1+.15)/C;high=(397-350-1-.05)/C
        us=np.clip(kp*(-.15+.02)+rp,0,rp)
        ue=np.clip(kp*(.05-.02)-rm,0,rp)
        self.assertEqual(us,0)
        self.assertGreater(C*k*low**2-1e-3,0)
        required=k*2.5*high/(cf-high)
        qmin=C*(cf-high)/.00055
        self.assertGreater(qmin*(ue-required)-1e-3,0)
        self.assertGreater(C*k*high**2-1e-3,0)

    def test_heat_selector_retains_strict_recovery_margin(self):
        C,k,NA,cb,V,cf,D,Qp,kap,h=sp.symbols('C k NA cb V cf D Qp kap h',positive=True)
        reaction=k*NA*cb;q=C*(cf-cb)/V
        b=D*k*NA/V*(cf-cb)-Qp
        a=-D*k*reaction*(NA+cb*V)/V
        ue=reaction/(cf-cb);uq=(-a+kap*h)/b
        nonnegative_numerator=q*Qp*ue+q*kap*h+Qp*C*k*cb**2
        self.assertEqual(sp.simplify(b*(q*(uq-ue)-C*k*cb**2)-nonnegative_numerator),0)


if __name__=='__main__':unittest.main()
