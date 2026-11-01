"""Selection/bookkeeping safety checks only: no event reads or classifier training."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
import zipfile

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import ggzh_rwt_manifest as M


def inventory(nominal, omitted=()):
    products=[]
    for sample,n in nominal.items():
        for seed in range(n['first_seed'],n['first_seed']+n['files']):
            if (sample,seed) in omitted: continue
            product=f'tables/{sample}-rwt10m-13600/table-{seed}.csv.gz'
            products.extend((product,product+'.sha256'))
    return dict(checked_utc='2026-11-01T00:00:00Z',products=products)


class Selection(unittest.TestCase):
    def test_all_available_files_are_kept_and_manifest_cannot_be_replaced(self):
        n=M.nominal()
        omitted=[('nlo-cp5',100872),('lo-cms',110001)]
        with tempfile.TemporaryDirectory() as td:
            directory=Path(td);i=directory/'inventory.json';o=directory/'manifest.json'
            i.write_text(json.dumps(inventory(n,omitted)))
            d=M.freeze(i,o)
            self.assertEqual(d['expected_events'],{'lo-cms':9995000,'nlo-cp5':9997500})
            self.assertEqual(len(list(M.products(d))),5998)
            self.assertEqual(d,M.load(o))
            with self.assertRaises(FileExistsError):M.freeze(i,o)
            d['nominal']['nlo-cp5']['files']-=1;d['selection_sha256']=M.digest(d)
            with self.assertRaises(ValueError):M.validate(d)

    def test_three_omissions_are_not_authorized(self):
        with tempfile.TemporaryDirectory() as td:
            directory=Path(td);i=directory/'inventory.json';o=directory/'manifest.json'
            i.write_text(json.dumps(inventory(M.nominal(),[('nlo-cp5',100871),('nlo-cp5',100872),('lo-cms',110001)])))
            with self.assertRaises(ValueError):M.freeze(i,o)
            self.assertFalse(o.exists())

    def test_loaded_files_event_counts_and_nominal_weights_must_match(self):
        with tempfile.TemporaryDirectory() as td:
            directory=Path(td);i=directory/'inventory.json';o=directory/'manifest.json'
            i.write_text(json.dumps(inventory(M.nominal(True))))
            d=M.freeze(i,o,True)
            for sample in M.SAMPLES:
                n=d['nominal'][sample]['events_per_file']
                meta=dict(files=[f'tables/{sample}/table-{s}.csv.gz' for s in d['included'][sample]],
                          events_per_file=[n]*2,events=n*2,nominal_column='2000' if sample=='nlo-cp5' else 'nominal')
                (directory/(sample+'.json')).write_text(json.dumps(meta))
            self.assertEqual(M.verify(d,directory),{'lo-cms':10000,'nlo-cp5':5000})
            p=directory/'nlo-cp5.json';meta=json.loads(p.read_text());meta['nominal_column']='nominal';p.write_text(json.dumps(meta))
            with self.assertRaises(ValueError):M.verify(d,directory)
            meta['nominal_column']='2000';p.write_text(json.dumps(meta))
            # A small NPY header fixture checks the cached array geometry without event data.
            header=repr(dict(descr='<f4',fortran_order=False,shape=(4999,))).encode()+b'\n'
            with zipfile.ZipFile(directory/'nlo-cp5.npz','w') as z:
                z.writestr('w0.npy',b'\x93NUMPY\x01\x00'+len(header).to_bytes(2,'little')+header)
            with self.assertRaisesRegex(ValueError,'array count differs'):M.verify(d,directory)


if __name__=='__main__':unittest.main()
