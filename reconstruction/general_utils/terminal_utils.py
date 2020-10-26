from subprocess import check_output, CalledProcessError
from tempfile import TemporaryFile


def get_out(*args):
    import getpass

    with TemporaryFile() as t:
        try:
            if args[0] == "chimerax":
              if(getpass.getuser() == "lcastillo"):
                args1 = list(args)
                args1[0] = "ChimeraX"
                args = tuple(args1)
            out = check_output(args, stderr=t)
            return 0, out.decode("utf-8")
        except CalledProcessError as e:
            t.seek(0)
            return e.returncode, t.read().decode("utf-8")
