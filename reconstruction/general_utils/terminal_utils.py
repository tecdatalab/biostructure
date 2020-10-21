from subprocess import check_output, CalledProcessError
from tempfile import TemporaryFile


def get_out(*args):
    with TemporaryFile() as t:
        try:
            out = check_output(args, stderr=t)
            return 0, out.decode("utf-8")
        except CalledProcessError as e:
            t.seek(0)
            return e.returncode, t.read().decode("utf-8")
