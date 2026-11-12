exec(open('tmp/storyboard_scenario_pilot/pilot.py').read().split('out={}')[0])
out={}
for q in (.15,1.,2.):
    out[f'trafficking_{q:g}']=run(f'trafficking_{q:g}',params={'ev_release':{'rab_conversion_baseline':.35*q,'escrt_baseline':.40*q,'rab_docking_baseline':.30*q}})
for tau in (.15,.5,1.5): out[f'reseal_{tau:g}']=run(f'reseal_{tau:g}',params={'ion_transport':{'tau_pore_reseal_s':tau}})
Path('tmp/storyboard_scenario_pilot/pathway_summary.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
