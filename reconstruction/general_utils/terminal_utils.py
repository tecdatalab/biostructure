def get_out(*args):
  import getpass, subprocess

  args1 = list(args)
  if args[0] == "chimerax":
    if (getpass.getuser() == "lcastillo"):
      args1[0] = "ChimeraX"

  result = subprocess.run(args1, stdout=subprocess.PIPE, check=True, timeout=900)
  exit_text = result.stdout.decode('utf-8')
  return result.returncode, exit_text
