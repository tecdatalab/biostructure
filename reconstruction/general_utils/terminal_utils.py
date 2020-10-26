

def get_out(*args):
    import getpass, subprocess

    args1 = list(args)
    if args[0] == "chimerax":
      if (getpass.getuser() == "lcastillo"):
        args1[0] = "ChimeraX"

    try:
      result = subprocess.run(args1, stdout=subprocess.PIPE)
      exit_text = result.stdout.decode('utf-8')
      return result.returncode ,exit_text

    except subprocess.CalledProcessError as grepexc:
      return grepexc.returncode, grepexc

