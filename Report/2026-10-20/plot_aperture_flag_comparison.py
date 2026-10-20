"""Plot fixed-denominator FINALCAT recovery before/after the flag ablation.

Fractions follow notebook 2.3 and Proposal: every unique FINALCAT row counts
once, overlap fields are unioned per band, and the denominator is the entire
valid-RM FINALCAT catalog. It is not restricted to the 7DT footprint. Existing
boundary masks are unchanged in the upstream catalog-production stage.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from astropy.io import fits
from astropy.table import Table

ROOT=Path(__file__).resolve().parents[2]
EXPERIMENT=ROOT/'output/experiments/aperture_flag_ablation'
BANDS=tuple(f'm{wave}' for wave in range(400,900,25))
FILENAME='master_catalog_MatchedWithFianlcat_OverallCosmosField.fits'


def selected_columns(path):
    """Restore only needed Boolean/numeric masks from the packed FITS masks."""
    with fits.open(path,memmap=True) as hdus:
        data=hdus['CATALOG'].data
        masks=hdus['COLUMN_MASKS'].data
        columns={name:i for i,name in enumerate(data.names)}
        bitorder=hdus['COLUMN_MASKS'].header['BITORDER']
        def value(name,fill=0):
            index=columns[name]
            bit=index%8 if bitorder=='little' else 7-index%8
            mask=((masks[:,index//8]>>bit)&1).astype(bool)
            result=np.array(data[name],copy=True)
            result[mask]=fill
            return result
        rm=value('finalcat_RM',np.nan)
        ids=value('finalcat_idx',-1)
        band_success=np.zeros((len(data),20),bool)
        for field in (1,2,3,4):
            prefix=f'7dt_cosmos{field}_'
            matched=value(prefix+'matched').astype(bool)
            for bi,band in enumerate(BANDS):
                band_success[:,bi]|=matched & value(prefix+band+'_success').astype(bool)
        counts=band_success.sum(1)
        assert np.array_equal(counts,value('n_success_bands'))
        assert len(np.unique(ids))==len(ids)
        return dict(rm=rm,ids=ids,success=band_success,nband=counts)


def main():
    output=EXPERIMENT/'comparison'
    output.mkdir(parents=True,exist_ok=True)
    samples={
        'baseline':selected_columns(EXPERIMENT/'baseline/output/singlephot'/FILENAME),
        # Once promoted, use the current branch's canonical combined product.
        'no_flag_cut':selected_columns((ROOT/'output/singlephot'/FILENAME)
            if (EXPERIMENT/'promotion.json').exists() else EXPERIMENT/'singlephot'/FILENAME),
    }
    reference=Table.read(ROOT/'COSMOS_data/COSMOS_finalcat.fits')
    reference_rm=np.asarray(np.ma.filled(reference['RM'],np.nan),float)
    valid=np.isfinite(reference_rm)&(reference_rm>0)&(reference_rm<50)
    # Use the original 0.5-mag bins and full-reference denominator, including
    # sources brighter than the displayed 16-21 mag interval.
    edges=np.arange(np.floor(reference_rm[valid].min()*2)/2,
                    np.ceil(reference_rm[valid].max()*2)/2+.5,.5)
    centers=(edges[:-1]+edges[1:])/2
    denominator=np.histogram(reference_rm[valid],bins=edges)[0]
    histograms={}
    cumulative={}
    for label,sample in samples.items():
        rm=sample['rm']; good=np.isfinite(rm)&(rm>0)&(rm<50)
        by_band=np.array([np.histogram(rm[good&sample['success'][:,bi]],bins=edges)[0] for bi in range(20)])
        any_counts=np.histogram(rm[good&(sample['nband']>0)],bins=edges)[0]
        histograms[label]=dict(band=by_band,any=any_counts)
        cumulative[label]={}
        for kind,selection in [('any',sample['nband']>0),('ge15',sample['nband']>=15),('all20',sample['nband']==20)]:
            cumulative[label][kind]=np.searchsorted(np.sort(rm[good&selection]),edges,side='right')
        for bi,band in enumerate(BANDS):
            cumulative[label][band]=np.searchsorted(np.sort(rm[good&sample['success'][:,bi]]),edges,side='right')
    def divide(numerator,den):
        return np.divide(numerator,den,out=np.full(np.shape(numerator),np.nan,float),where=den>0)
    wavelengths=np.arange(400,900,25)
    cmap=LinearSegmentedColormap.from_list('flag_ablation_wavelength',['#8A00CC','#003BFF','#00BFFF','#00A83B','#FFD400','#FF6500','#E00000'])
    norm=Normalize(400,875)
    fig,axes=plt.subplots(2,2,figsize=(14,10),sharex=True,sharey='row',layout='constrained',gridspec_kw={'height_ratios':[3,1]})
    for col,(label,title) in enumerate([('baseline','Baseline: aperture flag == 0'),('no_flag_cut','Aperture flag cut removed')]):
        ax,frac=axes[:,col]
        ax.plot(centers,denominator,'k--',lw=1.6,label='Full FINALCAT reference')
        for bi,band in enumerate(BANDS):
            color=cmap(norm(wavelengths[bi]))
            counts=histograms[label]['band'][bi]
            ax.plot(centers,counts,'o-',ms=2,lw=1,color=color)
            frac.plot(centers,divide(counts,denominator),'o-',ms=2,lw=1,color=color)
        ax.set(title=title,ylabel='Number of sources')
        ax.legend(fontsize=10)
        frac.set(xlabel='FINALCAT RM [AB mag]',ylabel='Band / reference',ylim=(0,1.05),xlim=(16,21))
        for a in (ax,frac):a.grid(alpha=.2)
    shown=(centers>=16)&(centers<=21)
    # Lines interpolate across the xlim boundaries. Include those boundary
    # heights so the rapidly rising reference curve is not clipped at RM=21.
    visible_reference=np.concatenate([denominator[shown],np.interp([16,21],centers,denominator)])
    axes[0,0].set_ylim(0,max(1,visible_reference.max())*1.10)
    scalar=plt.cm.ScalarMappable(norm=norm,cmap=cmap)
    fig.colorbar(scalar,ax=axes[0,:].tolist(),label='7DT central wavelength [nm]',fraction=.025)
    fig.suptitle('COSMOS1-4: identical raw photometry, footprint masks, reference and magnitude bins',fontsize=13)
    fraction_path=output/'fraction_before_after.png'
    fig.savefig(fraction_path,dpi=180)
    plt.close(fig)

    reference_cum=np.searchsorted(np.sort(reference_rm[valid]),edges,side='right')
    fig,axes=plt.subplots(1,2,figsize=(14,6),layout='constrained')
    for kind,color,title in [('any','tab:blue','Any successful band'),('ge15','tab:orange','At least 15 bands'),('all20','tab:green','All 20 bands')]:
        for label,ls in [('baseline','--'),('no_flag_cut','-')]:
            fraction=divide(cumulative[label][kind],reference_cum)
            axes[0].plot(edges,fraction,ls,color=color,label=f'{title}: {label}')
    axes[0].set(xlabel='RM upper limit [AB mag]',ylabel='Cumulative fraction of full FINALCAT',xlim=(16,21),ylim=(0,1.05))
    axes[0].legend(fontsize=8,loc='upper center',bbox_to_anchor=(.5,-.18),ncol=2)
    for label,ls in [('baseline','--'),('no_flag_cut','-')]:
        axes[1].plot(centers,divide(histograms[label]['any'],denominator),ls,color='tab:blue',label=label)
    axes[1].set(xlabel='FINALCAT RM [AB mag]',ylabel='Any-band / reference in bin',xlim=(16,21),ylim=(0,1.05))
    axes[1].legend()
    for ax in axes:ax.grid(alpha=.2)
    fig.suptitle('Dashed: baseline; solid: aperture flag cut removed')
    overall_path=output/'overall_fraction_before_after.png'
    fig.savefig(overall_path,dpi=180)
    plt.close(fig)

    rows=[]
    denominator21=int(np.sum(valid&(reference_rm<=21)))
    for bi,band in enumerate(BANDS):
        selections={}
        for label,sample in samples.items():
            use=np.isfinite(sample['rm'])&(sample['rm']>0)&(sample['rm']<=21)&sample['success'][:,bi]
            selections[label]=set(sample['ids'][use].tolist())
        before,after=selections['baseline'],selections['no_flag_cut']
        rows.append(dict(band=band,N_reference_RMle21=denominator21,N_before=len(before),N_after=len(after),
                         fraction_before=len(before)/denominator21,fraction_after=len(after)/denominator21,
                         increase_percentage_points=100*(len(after)-len(before))/denominator21,
                         gained_finalcat_rows=len(after-before),lost_finalcat_rows=len(before-after)))
    band_table=Table(rows=rows)
    band_table.write(output/'fractions_RMle21_by_band.csv',format='ascii.csv',overwrite=True)
    population_rows=[]
    for label,sample in samples.items():
        bright=np.isfinite(sample['rm'])&(sample['rm']>0)&(sample['rm']<=21)
        for name,use in [('any',sample['nband']>0),('ge15',sample['nband']>=15),('all20',sample['nband']==20)]:
            n=int(np.sum(bright&use))
            population_rows.append(dict(variant=label,population=name,RM_max=21,N_reference=denominator21,N_sources=n,fraction=n/denominator21))
    overall_table=Table(rows=population_rows)
    overall_table.write(output/'fractions_RMle21_populations.csv',format='ascii.csv',overwrite=True)
    bin_table=Table({'RM_lo':edges[:-1],'RM_hi':edges[1:],'N_reference':denominator})
    for label in samples:
        bin_table[label+'_N_any']=histograms[label]['any']
        bin_table[label+'_fraction_any']=divide(histograms[label]['any'],denominator)
        for bi,band in enumerate(BANDS):
            bin_table[label+'_N_'+band]=histograms[label]['band'][bi]
            bin_table[label+'_fraction_'+band]=divide(histograms[label]['band'][bi],denominator)
    bin_table.write(output/'fractions_by_RM_bin.csv',format='ascii.csv',overwrite=True)
    summary=dict(baseline_commit='d3c9272',branch='codex/aperture-flag-ablation',N_reference_RMle21=denominator21,
                 rows=[dict(zip(overall_table.colnames,[v.item() if isinstance(v,np.generic) else v for v in row])) for row in overall_table],
                 production_outputs_preserved=not (EXPERIMENT/'promotion.json').exists(),
                 historical_baseline_preserved=True,photometry_rerun=False,
                 denominator='All FINALCAT rows with finite positive RM <= threshold; no new footprint restriction.')
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(overall_table)
    print(band_table)
    return fraction_path,overall_path,overall_table,band_table


def save_proposal_outputs(results):
    """Persist the rerunnable comparison section and its actual rendered outputs.

    Proposal.ipynb is intentionally local-only under the repository policy;
    this versioned helper preserves the code needed to regenerate its section.
    """
    import base64
    import nbformat
    fraction_path,overall_path,populations,bands=results
    path=ROOT/'Notebooks/Proposal.ipynb'
    nb=nbformat.read(path,as_version=4)
    heading='## Aperture-flag ablation: fixed-denominator fraction comparison'
    prose=(heading+'\n\nBaseline commit: `d3c9272`; branch: `codex/aperture-flag-ablation`.\n'
           'Remove only `aper3_flags == 0`, then rerun COSMOS1–4 matching and downstream diagnostics. '
           'Reuse existing raw photometry; retain finite values, size and magnitude cuts, boundary masks, '
           'matching rules and photometric-success definitions.\n\n'
           'Baseline and experiment products are isolated under `output/experiments/aperture_flag_ablation/`. '
           'Both variants use the same full valid-RM FINALCAT denominator and 0.5-mag bins. '
           'Deduplicate overlapping fields by FINALCAT row and band. The tables use the exact inclusive '
           '`RM <= 21` threshold. Increased recovery alone does not establish photo-z accuracy.\n')
    code='''# Replot saved catalogs only; this cell never invokes SourceXtractor++.
from pathlib import Path
import importlib.util
from IPython.display import display, Image

ABLATION_ROOT = Path.cwd().resolve()
if ABLATION_ROOT.name == 'Notebooks':
    ABLATION_ROOT = ABLATION_ROOT.parent
ABLATION_HELPER = ABLATION_ROOT / 'Report/2026-10-20/plot_aperture_flag_comparison.py'
ABLATION_SPEC = importlib.util.spec_from_file_location('aperture_flag_comparison', ABLATION_HELPER)
ABLATION_MODULE = importlib.util.module_from_spec(ABLATION_SPEC)
ABLATION_SPEC.loader.exec_module(ABLATION_MODULE)
ABLATION_FRACTION_PNG, ABLATION_OVERALL_PNG, ABLATION_POPULATIONS, ABLATION_BANDS = ABLATION_MODULE.main()
display(Image(filename=str(ABLATION_FRACTION_PNG)))
display(Image(filename=str(ABLATION_OVERALL_PNG)))
display(ABLATION_POPULATIONS)
display(ABLATION_BANDS)
'''
    for i,cell in enumerate(nb.cells):
        if cell.cell_type=='markdown' and heading in cell.source:
            assert i+1<len(nb.cells) and nb.cells[i+1].cell_type=='code'
            nb.cells[i].source=prose
            target=nb.cells[i+1]
            target.source=code
            break
    else:
        nb.cells.append(nbformat.v4.new_markdown_cell(prose))
        target=nbformat.v4.new_code_cell(code)
        nb.cells.append(target)
    target.execution_count=max([c.get('execution_count') or 0 for c in nb.cells if c.cell_type=='code'])+1
    target.outputs=[]
    for image_path in (fraction_path,overall_path):
        target.outputs.append(nbformat.v4.new_output('display_data',data={
            'image/png':base64.b64encode(image_path.read_bytes()).decode('ascii'),
            'text/plain':f'<Image: {image_path.name}>'}))
    for table in (populations,bands):
        target.outputs.append(nbformat.v4.new_output('display_data',data={
            'text/html':table._repr_html_(),'text/plain':str(table)}))
    nbformat.validate(nb)
    nbformat.write(nb,path)
    print(f'Saved executed fraction section to {path}')


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--update-proposal',action='store_true')
    arguments=parser.parse_args()
    results=main()
    if arguments.update_proposal:
        save_proposal_outputs(results)
