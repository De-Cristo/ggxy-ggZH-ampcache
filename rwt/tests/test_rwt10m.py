"""Small schema and angular checks; no training and no event-sample reads."""
import gzip
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import ggzh_reweight_common as C
import ggzh_reweight_bdt as B


class Campaign(unittest.TestCase):
    def test_features_keep_all_closure_tests(self):
        for letter,features in C.FEATURE_SETS.items():
            B.configure(features)
            self.assertEqual(B.FEATURES,features)
            self.assertEqual(set(B.FEATURE_VARS),set(features))
            self.assertEqual(set(B.TEST_VARS),set(C.CLOSURE_VARS)-set(features))
            self.assertFalse({'pth','yh','drzh'} & set(features))
        self.assertEqual([len(C.FEATURE_SETS[k]) for k in 'ABC'],[6,9,17])

    def test_angles_wrap_and_missing_jet_stays_missing(self):
        cols={'phiz':np.array([np.pi-.1]),'phih':np.array([-np.pi+.1]),
              'ptz':np.array([50.]),'pth':np.array([50.]),'etaz':np.array([0.]),'etah':np.array([0.]),
              'ptj1':np.array([40.]),'etaj1':np.array([0.]),'phij1':np.array([np.pi]),
              'ptj2':np.array([-1.]),'etaj2':np.array([-1.]),'phij2':np.array([-1.])}
        C.append_derived(cols)
        self.assertAlmostEqual(float(cols['dphizzh'][0]),.1,places=6)
        self.assertAlmostEqual(float(cols['drj1zh'][0]),0.,places=6)
        self.assertEqual(float(cols['drj2zh'][0]),-1.)
        self.assertEqual(float(cols['dphij1j2'][0]),-1.)
        self.assertEqual(json.loads(C.dumps({'x':np.array([[np.nan,np.inf]])})),{'x':[[None,None]]})

    def test_table_schema_roundtrip_and_nominal(self):
        import ggzh_table as T
        base=['ievt','w:nominal']+['w:'+k for k in T.NANO_ORDER]+['w:2000']
        kin=['ptz','yz','mz','pth','yh','mzh','ptzh','yzh','drzh','etaz','phiz','etah','phih',
             'ptj1eta47','cell','njets30','njets30eta47']
        extra=['etaj1','phij1','ptj2','etaj2','phij2','ptl1','etal1','phil1','ptl2','etal2','phil2']
        with tempfile.TemporaryDirectory() as td:
            d=Path(td)
            for added in ([],extra):
                header=base+kin+added
                vals={k:1. for k in header};vals.update({'w:nominal':2.,'w:2000':3.,'mz':91.,'cell':0.})
                for k in T.NANO_ORDER: vals['w:'+k]=4.
                f=d/'table-100001.csv.gz'
                with gzip.open(f,'wt') as out:
                    out.write(','.join(header)+'\n')
                    for i in range(2):
                        vals['ievt']=i;out.write(','.join(str(vals[k]) for k in header)+'\n')
                path=d/'test.npz'
                T.main(['--label','test','--tables',str(f),'--out',str(path),'--nominal','2000'])
                with np.load(path) as out:
                    np.testing.assert_array_equal(out['w0'],[3.,3.])
                    np.testing.assert_array_equal(out['sw'],np.full((2,9),2.))
                    self.assertEqual('ptj2' in out, bool(added))

if __name__=='__main__': unittest.main()
