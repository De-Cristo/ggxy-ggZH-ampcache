#!/usr/bin/env python3
"""Launch a complete real ParT fold on an explicitly allocated external GPU worker.

No dependency installation/download or publication is performed. Input data are
a checksum-verified extraction of a frozen prepared bundle; source is snapshotted
before training. Each invocation creates a fresh work directory and preserves
failures, logs, runtime/environment information and model provenance.
"""
import argparse,hashlib,json,os,platform,re,resource,shutil,subprocess,sys,tarfile,time
from pathlib import Path

sys.dont_write_bytecode=True


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()


def unpack(archive,dest):
    with tarfile.open(archive) as t:
        assert all(not Path(m.name).is_absolute() and '..' not in Path(m.name).parts and not m.issym() and not m.islnk() for m in t.getmembers())
        t.extractall(dest)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--channel',choices=['zll','zvv'],required=True)
    ap.add_argument('--context',choices=['baseline24','analysis-jets-v1'],required=True)
    ap.add_argument('--task',type=int,required=True)
    ap.add_argument('--data',type=Path,required=True,help='Verified extracted parton-data directory')
    ap.add_argument('--work',type=Path,required=True,help='New independent directory; must not exist')
    ap.add_argument('--run-id',required=True,help='Unique external campaign ID; letters/digits/dot/underscore/hyphen')
    a=ap.parse_args()
    host=platform.node().split('.')[0]
    assert host!='cmscon' and not host.startswith('lxplus'),'Run only on an allocated compute worker'
    assert os.environ.get('RWT_EXTERNAL_COMPUTE_NODE')=='1','Declare the allocated compute worker explicitly'
    assert re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,127}',a.run_id)
    assert 0<=a.task<(5 if a.context=='analysis-jets-v1' else 15)
    assert os.environ.get('CUDA_VISIBLE_DEVICES','').strip(),'Select exactly one GPU through CUDA_VISIBLE_DEVICES'
    import torch
    assert torch.cuda.is_available() and torch.cuda.device_count()==1,'Exactly one assigned CUDA GPU required'
    source=Path(__file__).resolve().parents[1];data=a.data.resolve();work=a.work.resolve()
    assert data.is_dir() and (data/'SHA256SUMS').is_file()
    assert work!=source and source not in work.parents,'Put work outside the source checkout'
    # Prevent an inconsistent edit while git archive is produced. Dependencies are
    # already installed, so bytecode is disabled before importing code/launching.
    tracked=subprocess.run(['git','-C',str(source),'rev-parse','--show-toplevel'],text=True,capture_output=True)
    in_git=tracked.returncode==0 and Path(tracked.stdout.strip()).resolve()==source
    if in_git:
        assert not subprocess.check_output(['git','-C',str(source),'status','--porcelain'],text=True).strip(),'Commit code changes before training'
    else:assert re.fullmatch(r'[0-9a-f]{40}',(source/'VERSION').read_text().split()[0]),'Use the exported archive or a clean committed checkout'
    work.mkdir(parents=True,exist_ok=False);code=work/'code';code.mkdir()
    if in_git:
        archive=work/'source.tgz'
        subprocess.run(['git','-C',str(source),'archive','--format=tar.gz','-o',str(archive),'HEAD'],check=True)
        unpack(archive,code)
        expected_pin=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()
        assert (code/'VERSION').read_text().split()[0]==expected_pin,'VERSION must use git export-subst for this commit'
    else:
        # Copy only the code/manifest tree, never data or another run directory.
        for name in ('VERSION','.gitattributes','production','generator','shower','rwt','tools'):
            p=source/name
            if p.is_dir():shutil.copytree(p,code/name,ignore=shutil.ignore_patterns('__pycache__'))
            elif p.is_file():shutil.copy2(p,code/name)
    pin=(code/'VERSION').read_text().split()[0];assert re.fullmatch(r'[0-9a-f]{40}',pin)
    code_files=sorted(p for p in code.rglob('*') if p.is_file())
    (work/'code-SHA256SUMS').write_text(''.join(f'{digest(p)}  {p.relative_to(code)}\n' for p in code_files))
    # Every bulk content read occurs here, on the external compute worker.
    subprocess.run(['sha256sum','-c','--quiet','SHA256SUMS'],cwd=data,check=True)
    audit=json.loads((data/'audit.json').read_text())
    config_name='config-1m.json' if a.channel=='zll' else 'config-zvv-1m.json'
    cfg=json.loads((code/'rwt/parton'/config_name).read_text())
    assert audit['selection_sha256']==cfg['selection_sha256'] and audit['artifact_count']==3000
    assert all(v['rows']==1000000 and v['zero_overflow'] for v in audit['samples'].values())
    (code/'parton-data').symlink_to(data,target_is_directory=True)
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',RWT_EXECUTION_CONTEXT='external-gpu',RWT_EXTERNAL_RUN_ID=a.run_id,
        RWT_PARTON_CHANNEL=a.channel,RWT_PARTON_CONTEXT=a.context,RWT_THREADS='4',OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',MPLBACKEND='Agg')
    env.pop('_CONDOR_SCRATCH_DIR',None)
    (work/'environment.txt').write_text(subprocess.check_output([sys.executable,'-m','pip','freeze'],text=True))
    gpu=subprocess.run(['nvidia-smi','--query-gpu=index,name,memory.total,driver_version','--format=csv,noheader'],text=True,capture_output=True)
    (work/'gpu.txt').write_text(gpu.stdout+gpu.stderr)
    record=dict(run_id=a.run_id,channel=a.channel,context=a.context,task=a.task,package=pin,selection_sha256=cfg['selection_sha256'],
        execution_context='external-gpu',python=sys.version,torch=str(torch.__version__),visible_devices=env['CUDA_VISIBLE_DEVICES'],
        gpu=torch.cuda.get_device_name(0),data=str(data),code_file_hashes='code-SHA256SUMS',start_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
    (work/'external-run.json').write_text(json.dumps(record,indent=2)+'\n')
    start=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN)
    with (work/'train.log').open('w') as log:
        rc=subprocess.run(['/usr/bin/time','-v','-o',str(work/'resources.txt'),sys.executable,'-u','rwt/parton/fit.py',str(a.task)],cwd=code,env=env,stdout=log,stderr=subprocess.STDOUT).returncode
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    record.update(exit_code=rc,elapsed_seconds=time.monotonic()-start,cpu_seconds=(after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime),
                  allocated_gpu_wall_hours=(time.monotonic()-start)/3600)
    products=[]
    if rc==0:
        fit=json.loads((code/'parton-fit.json').read_text())
        assert fit['checkpoint_reload']=='PASS' and fit['selection_sha256']==cfg['selection_sha256']
        assert fit['task']==a.task and fit['external_run_id']==a.run_id and fit['package'].split()[0]==pin
        for name in ('parton-fit.tgz','parton-fit.json','parton-reload-audit.json'):
            p=code/name;assert p.stat().st_size>0
            p.with_name(p.name+'.sha256').write_text(digest(p)+'\n');products.append(name)
    record['verified_products']=products
    (work/'external-run.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2),flush=True)
    if rc:raise SystemExit(rc)

if __name__=='__main__':main()
