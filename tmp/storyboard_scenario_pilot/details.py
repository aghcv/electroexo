from electro_exocytosis.config import *
from electro_exocytosis.simulation import Simulation
import json
s=SimulationScenario(scenario=ScenarioConfig(name='control'),pulse=PulseConfig(amplitude_kV_cm=15.,pulse_width_ns=200.,pulse_number=20,repetition_rate_Hz=5.),exposure=ExposureConfig(dosimetry_model='joule_lumped_thermal'),simulation=SimulationConfig(t_end_s=1800.,output_dt_s=.25))
for nm,par in [('pump_low',{'ion_transport':{'SERCA_Vmax_uM_s':.0375,'PMCA_Vmax_uM_s':.05,'NCX_Vmax_uM_s':.015}}),('control',{}),('pump_high',{'ion_transport':{'SERCA_Vmax_uM_s':.3,'PMCA_Vmax_uM_s':.4,'NCX_Vmax_uM_s':.12}}),('reseal_fast',{'ion_transport':{'tau_pore_reseal_s':.15}}),('reseal_slow',{'ion_transport':{'tau_pore_reseal_s':1.5}})]:
 r=Simulation(s,params_override=par).run()
 print(nm,json.dumps({k:r.summary[k] for k in ['RNA_enrichment','protein_enrichment','lipid_enrichment','antigen_enrichment','potency_score','cumulative_small_EV','cumulative_medium_large_EV','cumulative_apoptotic_body','viability_fraction','bona_fide_EV_fraction']}))
for q in (.25,1.,2.):
 r=Simulation(s,params_override={'ev_release':{'rab_conversion_baseline':.35*q,'escrt_baseline':.40*q,'rab_docking_baseline':.30*q}}).run()
 print('routing_'+str(q),json.dumps({k:r.summary[k] for k in ['RNA_enrichment','protein_enrichment','lipid_enrichment','antigen_enrichment','potency_score','peak_secretory_bias','cumulative_small_EV','cumulative_medium_large_EV','cumulative_apoptotic_body','viability_fraction','bona_fide_EV_fraction']}))
 print('terminal', {k:r.parameters_used['terminal_ev_release'][k] for k in ['escrt_dependent_signal','rab_conversion_signal','rab_docking_signal','secretory_bias']})
