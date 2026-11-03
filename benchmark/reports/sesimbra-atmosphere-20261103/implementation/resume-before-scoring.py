"""Resume exports/glint after wrapper fixes; atmospheric arrays are untouched."""
from process import *
for branch in ['s2_radcor_fixed','s2_radcor_auto','pneo_marine_control','pneo_continental']:
 scene=S2 if branch.startswith('s2') else PN
 out=ROOT/branch
 result=next(out.glob('*L2R.nc'))
 wanted=(442,492,559,665,833,1610) if scene==S2 else (439,486,563,655,826)
 ds=[b+'_'+str(w) for b in ('rhos','rhosu') for w in wanted]
 export(result,scene,ds)
 if scene==S2:
  glint=out/'glint'
  complete=list(glint.glob('*rhos_*.tif')) if glint.exists() else []
  if not complete:
   if glint.exists():glint.rename(out/'glint-initial-aborted')
   glint.mkdir()
   settings=dict(ac.settings['defaults'])
   settings.update(ac.acolite.settings.parse('S2B_MSI'))
   settings.update(output=str(glint),glint_write_rhog_ref=True,dsf_interface_reflectance=False,ancillary_data=False)
   ac.settings['run']=settings.copy();ac.settings['user']=settings.copy()
   ac.glint.default(str(result),settings=settings,new_file=True)
  gf=list(glint.glob('*.nc'));assert len(gf)==1
  export(gf[0],scene,ds)
 with Dataset(result) as nc:
  gatts={k:nc.getncattr(k) for k in nc.ncattrs() if k in ['sensor','sza','vza','raa','ac_aot_550','ac_model','radcor_version']}
 dump(out/'completed.json',{'branch':branch,'l2r':str(result),'metadata':gatts,'grid_matches_source':True,
   'l1r_sha256':sha(ROOT/'inputs'/scene/(scene+'_L1R.nc')),'completion_source_sha256':sha(Path(__file__))})
 print('Finalized',branch,gatts,flush=True)
