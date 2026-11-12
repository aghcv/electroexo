from electro_exocytosis.config import *
from electro_exocytosis.simulation import Simulation
from examples.storyboard_multilayer_overview import SCENARIOS as REGIMES
import time, json
from pathlib import Path

basepulse=dict(amplitude_kV_cm=15., pulse_width_ns=200., pulse_number=20, repetition_rate_Hz=5.)
def run(name,pulse=None,cell=None,params=None):
    s=SimulationScenario(scenario=ScenarioConfig(name=name),pulse=PulseConfig(**(basepulse|(pulse or {}))), exposure=ExposureConfig(dosimetry_model='joule_lumped_thermal'),cell_state=CellStateConfig(**(cell or {})),simulation=SimulationConfig(t_end_s=1800.,output_dt_s=.25))
    t=time.time(); r=Simulation(s,params_override=params or {}).run(); print(name,round(time.time()-t,2),flush=True)
    d={k:r.summary[k] for k in ('peak_ca_i','peak_ros','min_atp','cumulative_small_EV','cumulative_medium_large_EV','cumulative_apoptotic_body','viability_fraction','bona_fide_EV_fraction','purity_score','potency_score','batch_adjusted_yield','peak_repair_state')}
    d.update(deltaVm=r.parameters_used['electrodynamics']['delta_Vm'],permeability=r.parameters_used['electrodynamics']['membrane_permeability'],ion_duration=r.parameters_used['effective_ion_perturbation_s'],peak_temperature=37+r.parameters_used['dosimetry']['temperature_rise_K'])
    d.update(time_ca_above_1=float(.25*(r.state_timeseries['Ca_i']>1).sum()),peak_damage=float(r.state_timeseries['damage'].max()))
    r.state_timeseries.to_csv(f'tmp/storyboard_scenario_pilot/{name}_state.csv',index=False)
    return d

out={}
for sp in REGIMES: out[sp.name]=run(sp.name,sp.pulse_overrides,sp.cell_state_overrides,sp.parameter_overrides)
for n in (5,20,80): out[f'number_{n:g}']=run(f'number_{n:g}',{'pulse_number':n})
for q in (.4,1.,2.5): out[f'cell_ca_{q:g}']=run(f'cell_ca_{q:g}',cell={'calcium_handling_modifier':q})
for q in (.25,1.,2.): out[f'pump_{q:g}']=run(f'pump_{q:g}',params={'ion_transport':{'SERCA_Vmax_uM_s':.15*q,'PMCA_Vmax_uM_s':.2*q,'NCX_Vmax_uM_s':.06*q}})
for q in (.25,1.,2.): out[f'repair_{q:g}']=run(f'repair_{q:g}',params={'remodeling_repair':{'resealing_annexin_weight':.50*q,'resealing_lysosome_weight':.35*q}})
for q in (.1,.3,.6): out[f'recovery_{q:g}']=run(f'recovery_{q:g}',params={'manufacturing_qc':{'isolation_efficiency':q}})
Path('tmp/storyboard_scenario_pilot/summary.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
