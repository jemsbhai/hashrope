"""Provenance/integrity contracts written before the fixture generator."""
import hashlib,json,pathlib,subprocess,sys,tempfile,unittest
ROOT=pathlib.Path(__file__).resolve().parents[1]
WHEEL=ROOT/'oracle'/'hashrope-0.2.2-py3-none-any.whl'
PIN='710112e3aa14a991918a15331c1f418ce9d1c9d93c34a9da09e22be33c68a2d2'
class GeneratorContract(unittest.TestCase):
 def test_pinned_wheel_origin_and_corpus_integrity(self):
  with tempfile.TemporaryDirectory() as tmp:
   out=pathlib.Path(tmp)/'fixtures'
   p=subprocess.run([sys.executable,'-I',str(ROOT/'scripts'/'generate-oracle.py'),'--output',str(out),'--wheel',str(WHEEL)],capture_output=True,text=True)
   self.assertEqual(p.returncode,0,p.stderr)
   m=json.loads((out/'manifest.json').read_text())
   data=(out/'oracle.json').read_bytes(); corpus=json.loads(data)
   self.assertEqual(m['wheel']['sha256'],PIN)
   self.assertTrue(all('.whl/' in p.replace('\\','/') for p in m['importOrigins'].values()))
   self.assertEqual(m['corpus']['sha256'],hashlib.sha256(data).hexdigest())
   self.assertEqual(m['corpus']['records'],len(corpus['rows']))
   self.assertEqual(len({r['id'] for r in corpus['rows']}),len(corpus['rows']))
   self.assertGreater(len(corpus['rows']),3000)
   self.assertEqual(len(m['pythonExports']),21)
   p2=subprocess.run([sys.executable,'-I',str(ROOT/'scripts'/'generate-oracle.py'),'--output',str(out),'--wheel',str(WHEEL)],capture_output=True,text=True)
   self.assertEqual(p2.returncode,0,p2.stderr)
   self.assertEqual(data,(out/'oracle.json').read_bytes())
 def test_corrupt_wheel_is_rejected_before_import(self):
  with tempfile.TemporaryDirectory() as tmp:
   bad=pathlib.Path(tmp)/'hashrope-0.2.2-py3-none-any.whl'; bad.write_bytes(WHEEL.read_bytes()+b'corrupt')
   p=subprocess.run([sys.executable,'-I',str(ROOT/'scripts'/'generate-oracle.py'),'--output',str(pathlib.Path(tmp)/'out'),'--wheel',str(bad)],capture_output=True,text=True)
   self.assertNotEqual(p.returncode,0)
   self.assertIn('wheel SHA256 mismatch',p.stderr)
   self.assertFalse((pathlib.Path(tmp)/'out'/'oracle.json').exists())
if __name__=='__main__': unittest.main()
