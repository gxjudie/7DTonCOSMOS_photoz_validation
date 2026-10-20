"""Promote no-aperture-flag products and persist downstream notebook outputs.

Raw SE++ catalogs and boundary masks are immutable inputs. Per-field executed
notebooks are stored with the experiment; canonical notebooks show COSMOS1
and the updated cross-field results. Each cell is saved immediately so a
long-running depth calculation never loses earlier outputs.
"""
from pathlib import Path
import argparse
import json
import os
import shutil
import time
import nbformat
from nbclient import NotebookClient

ROOT=Path(__file__).resolve().parents[2]
EXP=ROOT/'output/experiments/aperture_flag_ablation'
EXECUTED=EXP/'executed_notebooks'
NAMES=['2_singphotandmatching.ipynb','2.1_singlephotcatalog_visualization.ipynb',
       '2.2_7DTdepthAnalysis.ipynb','2.3_numofmatchedbands_digonostics.ipynb','Proposal.ipynb']


def promote():
    """Archive all old derived files before replacing canonical product trees."""
    archive=EXP/'baseline/full_products'
    assert not archive.exists(), 'Promotion already started; do not overwrite the baseline.'
    archive.mkdir(parents=True)
    for name in NAMES:
        target=EXP/'pre_refresh_notebooks'/name
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(ROOT/'Notebooks'/name,target)
    for field in (1,2,3,4):
        canonical=ROOT/f'output/COSMOS{field}/singlephot'
        snapshot=archive/f'COSMOS{field}/singlephot'
        snapshot.mkdir(parents=True)
        for dirname in ('diagnostics','analysis'):
            old=canonical/dirname
            if old.exists():
                # Move, rather than copy, to avoid carrying stale source-ID
                # figures into the new branch products and save disk space.
                shutil.move(str(old),str(snapshot/dirname))
            shutil.copytree(EXP/f'COSMOS{field}/singlephot'/dirname,old)
        master=f'cosmos{field}_master_catalog.fits'
        shutil.copy2(canonical/master,snapshot/master)
        shutil.copy2(EXP/f'COSMOS{field}/singlephot'/master,canonical/master)
    canonical=ROOT/'output/singlephot'
    shutil.copytree(canonical,archive/'singlephot')
    # Retain no catalog-dependent files from the baseline combined tree.
    shutil.rmtree(canonical)
    shutil.copytree(EXP/'singlephot',canonical)
    EXECUTED.mkdir(parents=True,exist_ok=True)
    (EXP/'promotion.json').write_text(json.dumps(dict(
        condition='aperture flags unrestricted; all other source cuts unchanged',
        baseline_archive=str(archive),canonical_products_updated=True,
        promoted_from=str(EXP)),indent=2)+'\n')
    print('Promoted all four fields and the combined catalog.',flush=True)


def execute(name,field=None,per_field_only=False,canonical=True):
    """Run a real photoz Jupyter kernel and save rich outputs after every cell."""
    os.chdir(ROOT)
    original=ROOT/'Notebooks'/name
    nb=nbformat.read(original,as_version=4)
    # Cell outputs from a different source-ID assignment must never survive.
    for cell in nb.cells:
        if cell.cell_type=='code':
            cell.outputs=[]
            cell.execution_count=None
    if field is not None:
        for cell in nb.cells:
            if cell.cell_type=='code':
                cell.source=cell.source.replace('COSMOS = 1\n',f'COSMOS = {field}\n',1)
    if per_field_only:
        # The first section contains all field-specific depth and magnitude
        # diagnostics. Combined sky/recovery sections run once in the primary.
        nb.cells=nb.cells[:23]
    nb.cells.insert(0,nbformat.v4.new_markdown_cell(
        '# Executed branch result: aperture flags unrestricted\n\n'
        'Existing SE++ photometry; original boundary masks and remaining cuts. '
        'Catalogs under the standard `output/` paths use the branch condition. '
        'Historical flag-zero results remain in the experiment baseline backup.'))
    nb.cells.insert(1,nbformat.v4.new_code_cell(
        "# Capture all scientific figures in this saved notebook.\n"
        "get_ipython().run_line_magic('matplotlib', 'inline')\n"
        "from pathlib import Path\nimport subprocess\n"
        "_original_subprocess_run = subprocess.run\n"
        "def _reject_source_extraction(*args, **kwargs):\n"
        "    command = args[0] if args else kwargs.get('args', '')\n"
        "    if 'sourcextractor' in str(command).lower():\n"
        "        raise RuntimeError('SE++ execution disabled: reuse existing raw catalogs.')\n"
        "    return _original_subprocess_run(*args, **kwargs)\n"
        "subprocess.run = _reject_source_extraction\n"))
    tag=f'COSMOS{field}_' if field else ''
    destination=EXECUTED/(tag+name)
    client=NotebookClient(nb,timeout=10800,kernel_name='photoz',
                          resources={'metadata':{'path':str(ROOT)}},
                          allow_errors=False,record_timing=True)
    with client.setup_kernel():
        for index,cell in enumerate(nb.cells):
            if cell.cell_type!='code':
                continue
            preview=cell.source.splitlines()[0][:95] if cell.source else ''
            print(f'{tag}{name} CELL {index}/{len(nb.cells)-1} START {time.strftime("%H:%M:%S")} {preview}',flush=True)
            try:
                client.execute_cell(cell,index)
            except BaseException:
                nbformat.write(nb,destination)
                print(f'FAILED: saved partial notebook to {destination}',flush=True)
                raise
            nbformat.write(nb,destination)
            print(f'CELL {index} DONE outputs={len(cell.outputs)}',flush=True)
    if canonical:
        # Keep original cell indexing/source layout in canonical notebooks.
        # Only the rich results and intended source edits should remain.
        nb.cells=nb.cells[2:]
        nbformat.write(nb,original)
        if name == 'Proposal.ipynb':
            save_versioned_proposal(nb)
    print(f'COMPLETE {destination}',flush=True)



def save_versioned_proposal(notebook):
    """Preserve the local-only Proposal result in a versioned numbered copy."""
    import copy
    snapshot=copy.deepcopy(notebook)
    snapshot.metadata['7dt_selection']={
        'branch':'codex/aperture-flag-ablation','aperture_flag_cut':False,
        'raw_photometry_reused':True,'boundary_masks_unchanged':True,
        'fields':[1,2,3,4],'execution_date':'2026-10-20'}
    note='> Executed branch snapshot: aperture flags unrestricted; original boundary masks and remaining selection criteria retained. Existing SE++ measurements are reused.'
    if note not in snapshot.cells[0].source:
        snapshot.cells[0].source+='\n\n'+note+'\n'
    nbformat.write(snapshot,ROOT/'Notebooks/4_Proposal_aperture_flag_results.ipynb')


def rerun_cells(name, indices):
    """Refresh selected independent sections while preserving other rich outputs."""
    os.chdir(ROOT)
    path = ROOT / 'Notebooks' / name
    notebook = nbformat.read(path, as_version=4)
    context = nbformat.v4.new_code_cell(
        "get_ipython().run_line_magic('matplotlib', 'inline')\n")
    notebook.cells.insert(0, context)
    destination = EXECUTED / ('COSMOS1_' + name)
    client = NotebookClient(notebook, timeout=10800, kernel_name='photoz',
                            resources={'metadata': {'path': str(ROOT)}},
                            allow_errors=False, record_timing=True)
    with client.setup_kernel():
        client.execute_cell(context, 0)
        for index in indices:
            print(f'{name} REFRESH CELL {index} START {time.strftime("%H:%M:%S")}', flush=True)
            cell = notebook.cells[index + 1]
            try:
                client.execute_cell(cell, index + 1)
            finally:
                saved = nbformat.from_dict(dict(notebook))
                saved.cells = notebook.cells[1:]
                nbformat.write(saved, path)
                nbformat.write(saved, destination)
            print(f'CELL {index} DONE outputs={len(cell.outputs)}', flush=True)
    print(f'COMPLETE {destination}', flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['promote','field_catalogs','field_visuals','field_depths','primary_depth','finalize_depth','combined','proposal'])
    args=parser.parse_args()
    if args.stage=='promote':
        promote()
    elif args.stage=='field_catalogs':
        for field in (1,2,3,4):
            execute(NAMES[0],field,canonical=field==1)
    elif args.stage=='field_visuals':
        for field in (1,2,3,4):
            execute(NAMES[1],field,canonical=field==1)
    elif args.stage=='field_depths':
        for field in (2,3,4):
            execute(NAMES[2],field,per_field_only=True,canonical=False)
    elif args.stage=='primary_depth':
        execute(NAMES[2],1)
    elif args.stage=='finalize_depth':
        rerun_cells(NAMES[2],[2,4,10,24,25,45,47,48])
    elif args.stage=='combined':
        execute(NAMES[3])
    else:
        execute(NAMES[4])


if __name__=='__main__':
    main()
