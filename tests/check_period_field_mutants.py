"""Compile deliberate field regressions with ICK in temporary source copies.

Run on a clean, approved candidate checkout with the pinned ICK stage.
Compilation failure is NOT a mutant kill: each mutant must build and reach
an assertion in the period test. No repository, branch or compiler is edited.
"""
import resource
import shutil
import signal
import subprocess
import sys
import tempfile
from pathlib import Path


MUTANTS = (
    ('wrong-work-field', 'src/periods.c',
     'exact_work work = {.field = surface->polygon.field};',
     'exact_work work = {.field = BILLIARD_SQRT3_TWELFTHS};',
     b'billiard_periods(&eighth,&eighth_ledger)'),
    ('missing-ledger-identity-check', 'src/periods.c',
     'ledger->field != rebuilt.field || ', '',
     b'billiard_periods_validate(surface,&changed)'),
    ('wrong-serialized-radicand', 'src/periods.c',
     'ledger->field == BILLIARD_SQRT2_EIGHTHS ? 2 : 3',
     'ledger->field == BILLIARD_SQRT2_EIGHTHS ? 3 : 2', b'strcmp(line'),
    ('wrong-affine-display-field', 'src/affine.c',
     'billiard_polygon_display(field,point.x)',
     'billiard_exact_display(point.x)', b'Assertion'),
    ('wrong-torus-determinant-field', 'src/affine.c',
     'exact_work work = {.field = periods->field};\n    exact_number determinant',
     'exact_work work = {.field = BILLIARD_SQRT3_TWELFTHS};\n    exact_number determinant',
     b'billiard_affine_apply(surface,&ledger,actions[action],&deformed)'),
)


def main():
    if len(sys.argv) != 3:
        raise SystemExit('usage: check_period_field_mutants.py ICK_CC ICK_SOURCE_DIR')
    root = Path(__file__).resolve().parent.parent
    cc, source = (str(Path(arg).absolute()) for arg in sys.argv[1:])
    revision = subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'], text=True).strip()
    # Assertion failures must not leave potentially large core files.
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    for name, path, old, new, expected in (('baseline',None,None,None,None), *MUTANTS):
        with tempfile.TemporaryDirectory(prefix='billiard-period-') as temp:
            dest = Path(temp)
            for directory in ('src','tests'):
                shutil.copytree(root/directory,dest/directory)
            shutil.copy2(root/'Makefile',dest/'Makefile')
            if path is not None:
                target = dest/path
                original = target.read_text()
                if original.count(old)!=1:
                    raise RuntimeError(f'{name}: source mutation anchor not unique')
                target.write_text(original.replace(old,new,1))
            # The existing Makefile check-ick enforces the actual compiler pins.
            subprocess.run(['make','-C',str(dest),'build/test_periods',
                            f'ICK_CC={cc}',f'ICK_SOURCE_DIR={source}',
                            f'SOURCE_REVISION={revision}','SOURCE_STATE=working_tree'],check=True)
            exe = dest/'build/test_periods'
            result = subprocess.run([str(exe)],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            if path is None:
                if result.returncode:
                    raise RuntimeError(f'baseline failed: {result.stderr!r}')
                subprocess.run([sys.executable,str(dest/'tests/check_eighth_periods.py'),str(exe)],check=True)
                print('baseline: PASS',flush=True)
            else:
                if result.returncode != -signal.SIGABRT or b'Assertion' not in result.stderr or expected not in result.stderr:
                    raise RuntimeError(f'{name}: not killed by expected assertion: {result.returncode}, {result.stderr!r}')
                print(f'{name}: rejected by runtime assertion',flush=True)
    print('5 compiled field mutants rejected')


if __name__=='__main__':
    main()
