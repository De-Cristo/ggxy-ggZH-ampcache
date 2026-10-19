"""The local executor must carry the reweighting code into each job directory.

Uses a tiny fake production job; no generator, shower or training is run.
Run with: python3 -m unittest discover -s production/tests
"""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


class LocalReweightingJob(unittest.TestCase):
    def test_reweighting_code_reaches_job(self):
        for archived in (True, False):
            with self.subTest(git_archive=archived), tempfile.TemporaryDirectory() as td:
                root = Path(td)
                pkg = root / "package"
                for name in ("production", "generator", "shower", "rwt"):
                    (pkg / name).mkdir(parents=True)
                    (pkg / name / ".keep").touch()
                (pkg / "VERSION").write_text("test\n")
                shutil.copy2(Path(__file__).parents[1] / "run_local.sh", pkg / "production/run_local.sh")
                (pkg / "rwt/run.sh").write_text("#!/bin/bash\nprintf 'model produced\\n' > model.txt\n")
                (pkg / "production/paper.sh").write_text('''#!/bin/bash
case $1 in
  list) echo 'rwt-model 1 1 1 1' ;;
  inputs) ;;
  outputs) echo 'model.txt rwt/model-vv-13600.txt' ;;
  run) bash rwt/run.sh ;;
  *) exit 2 ;;
esac
''')
                (pkg / "production/paper.sh").chmod(0o755)
                env = os.environ.copy()
                if archived:
                    def git(*args):
                        subprocess.run(["git", *args], cwd=pkg, check=True, stdout=subprocess.DEVNULL)
                    git("init", "-q")
                    git("add", ".")
                    git("-c", "user.name=Test", "-c", "user.email=test@example.invalid",
                        "commit", "-qm", "fixture")
                else:
                    # Force the source-distribution fallback even if the test's parent is a repository.
                    (root / "bin").mkdir()
                    fakegit = root / "bin/git"
                    fakegit.write_text("#!/bin/bash\nexit 1\n")
                    fakegit.chmod(0o755)
                    env["PATH"] = str(root / "bin") + os.pathsep + env["PATH"]
                subprocess.run(["bash", str(pkg / "production/run_local.sh"), "vv13600", str(root / "work")],
                               check=True, env=env, stdout=subprocess.DEVNULL)
                self.assertEqual((root / "work/products/rwt/model-vv-13600.txt").read_text(), "model produced\n")
                self.assertTrue((root / "work/jobs/vv13600/rwt-model/0/.done").is_file())


if __name__ == "__main__":
    unittest.main()
