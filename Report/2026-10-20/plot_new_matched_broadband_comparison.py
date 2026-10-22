"""Compare new matched 7DT detections to FINALCAT GM/RM in COSMOS1-4.

Pool field measurements without deduplicating overlapping sky objects, as
requested. GM and RM are FINALCAT magnitudes, never SDSS g/r. Retain every
valid saved band measurement, without S/N, color, aperture-flag or outlier
cuts. The medium/broad-band comparison is descriptive, not a zero-point fit.
"""
from pathlib import Path
import numpy as np
from astropy.table import Table
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from scl_plot import scl_plt

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'output/singlephot/analysis'
STEM = 'new_matched_sources_m475_vs_g_m600_vs_r'


def values(table, name, fill=np.nan):
    """Honor native FITS masks, including absent-band integer source IDs."""
    return np.asarray(np.ma.filled(table[name],fill))


def main():
    rows, count_rows = [], []
    for field in (1,2,3,4):
        path = ROOT/'output'/f'COSMOS{field}_flagcut_removed_new_7DT_sources_matched_to_FINALCAT.fits'
        catalog = Table.read(path)
        assert np.all(catalog['is_entirely_new_7dt_source'])
        for band, reference in [('m475','GM'),('m600','RM')]:
            if band+'_mag' not in catalog.colnames:
                count_rows.append(dict(COSMOS=field,band=band,N=0))
                continue
            x = values(catalog,'finalcat_'+reference)
            y = values(catalog,band+'_mag')
            error = values(catalog,band+'_mag_err')
            source_id = values(catalog,band+'_source_id',-1)
            # Reject missing/sentinel values only. A positive finite mag error
            # is required for a valid saved detection; its size is unrestricted.
            good = (np.isfinite(x)&(x>0)&(x<50)&np.isfinite(y)&(y>0)&(y<50)
                    &np.isfinite(error)&(error>0)&(source_id>=0))
            count_rows.append(dict(COSMOS=field,band=band,N=int(good.sum())))
            for i in np.flatnonzero(good):
                rows.append(dict(COSMOS=field,master_id=int(catalog['master_id'][i]),
                    finalcat_idx=int(catalog['finalcat_idx'][i]),band=band,
                    finalcat_column=reference,finalcat_mag=float(x[i]),
                    mag_7dt=float(y[i]),mag_err_7dt=float(error[i]),
                    source_id=int(source_id[i]),aperture_flag=int(catalog[band+'_aper3_flags'][i])))
    points = Table(rows=rows)
    counts = Table(rows=count_rows)
    OUT.mkdir(parents=True,exist_ok=True)
    points.write(OUT/(STEM+'_points.csv'),format='ascii.csv',overwrite=True)
    with plt.rc_context():
        scl_plt.default()
        # These sizes remain legible when the 10.5-inch figure is reduced to
        # a typical two-column paper width. Embed TrueType fonts in the PDF.
        plt.rcParams.update({'axes.labelsize':24,'xtick.labelsize':19,
                             'ytick.labelsize':19,'pdf.fonttype':42,'ps.fonttype':42})
        fig, axes = plt.subplots(1,2,figsize=(10.5,5.0),layout='constrained')
        for ax,band,letter in zip(axes,('m475','m600'),('g','r')):
            selected = points[np.asarray(points['band']).astype(str)==band]
            x,y = np.asarray(selected['finalcat_mag']),np.asarray(selected['mag_7dt'])
            # Include every valid point on both axes; no percentile clipping.
            lo = np.floor((min(x.min(),y.min())-.2)*2)/2
            hi = np.ceil((max(x.max(),y.max())+.2)*2)/2
            ax.plot([lo,hi],[lo,hi],color='0.35',ls='--',lw=1.4,zorder=1)
            ax.scatter(x,y,s=22,color='#2878B5',alpha=.65,linewidths=0,zorder=2)
            ax.set(xlabel=f'${letter}$',ylabel=rf'$m_{{{band[1:]}}}$',
                   xlim=(lo,hi),ylim=(lo,hi))
            ax.set_aspect('equal',adjustable='box')
            ax.xaxis.set_major_locator(MultipleLocator(2))
            ax.yaxis.set_major_locator(MultipleLocator(2))
            ax.xaxis.set_minor_locator(MultipleLocator(.5))
            ax.yaxis.set_minor_locator(MultipleLocator(.5))
            assert ax.get_title()==''
        pdf, png = OUT/(STEM+'.pdf'),OUT/(STEM+'.png')
        fig.savefig(pdf)
        fig.savefig(png,dpi=450)
        plt.close(fig)
    print(counts)
    for band in ('m475','m600'):
        sample=points[np.asarray(points['band']).astype(str)==band]
        print(band,'field measurements',len(sample),'unique FINALCAT rows',len(np.unique(sample['finalcat_idx'])))
    print(pdf);print(png)
    return pdf,png,counts


def save_notebook_outputs(results):
    """Keep the rerunnable plotting cell and its actual figure in notebooks."""
    import base64
    import nbformat
    pdf,png,counts=results
    heading='## New matched 7DT detections versus FINALCAT g and r'
    prose=heading+'''\n\nPool COSMOS1-4 field measurements from the strictly new 7DT + FINALCAT\nmatched subset. The horizontal coordinates are FINALCAT **GM and RM**,\nnot SDSS magnitudes. Require finite 0<mag<50, a saved source ID, and a finite\npositive 7DT magnitude error. Apply no S/N, flag or outlier cut; retain overlap\nfield measurements separately. COSMOS4 has no m475 image. Dashed lines show\nidentity, not a fitted relation. Medium- and broad-band filters differ, so\ndeviations alone do not identify a zero-point error. Use scl_plot defaults,\nlarge labels, no titles, vector PDF and 450-dpi PNG.\n'''
    code='''# Recreate the figure using saved catalogs only.
from pathlib import Path
import importlib.util
from IPython.display import display, Image
BB_ROOT = Path.cwd().resolve()
if BB_ROOT.name == 'Notebooks':
    BB_ROOT = BB_ROOT.parent
BB_SPEC = importlib.util.spec_from_file_location(
    'new_matched_broadband_plot',
    BB_ROOT / 'Report/2026-10-20/plot_new_matched_broadband_comparison.py')
BB_MODULE = importlib.util.module_from_spec(BB_SPEC)
BB_SPEC.loader.exec_module(BB_MODULE)
BB_PDF, BB_PNG, BB_COUNTS = BB_MODULE.main()
display(Image(filename=str(BB_PNG)))
display(BB_COUNTS)
'''
    cell=nbformat.v4.new_code_cell(code)
    cell.outputs=[nbformat.v4.new_output('display_data',data={
        'image/png':base64.b64encode(png.read_bytes()).decode(),
        'text/plain':'<New matched sources: m475 versus g, m600 versus r>'}),
        nbformat.v4.new_output('display_data',data={'text/html':counts._repr_html_(),'text/plain':str(counts)})]
    for name in ['Proposal.ipynb','4_Proposal_aperture_flag_results.ipynb','2.3_numofmatchedbands_digonostics.ipynb']:
        path=ROOT/'Notebooks'/name
        notebook=nbformat.read(path,4)
        for i,c in enumerate(notebook.cells):
            if c.cell_type=='markdown' and c.source.startswith(heading):
                del notebook.cells[i:i+2]
                break
        notebook.cells.extend([nbformat.v4.new_markdown_cell(prose),cell])
        nbformat.validate(notebook)
        nbformat.write(notebook,path)
        print('Saved broadband comparison:',path)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--update-notebooks',action='store_true')
    args=parser.parse_args()
    result=main()
    if args.update_notebooks:
        save_notebook_outputs(result)
