import subprocess


def get_out(*args):
  import getpass, subprocess

  args1 = list(args)
  if args[0] == "chimerax":
    if (getpass.getuser() == "lcastillo"):
      args1[0] = "ChimeraX"

  result = subprocess.run(args1, stdout=subprocess.PIPE, check=True, timeout=900)
  exit_text = result.stdout.decode('utf-8')
  return result.returncode, exit_text

  # cmd = ""
  # for i in args:
  #   cmd += i + " "
  #
  # execute_command(cmd)


# Excecute command
def execute_command(cmd, printExit=False):
  try:
    if printExit:
      subprocess.check_call([cmd], shell=True)
    else:
      subprocess.check_call([cmd],
                            shell=True,
                            stdout=subprocess.DEVNULL,
                            stderr=subprocess.STDOUT)

  except subprocess.CalledProcessError:
    raise RuntimeError('Command "%s" does not work' % cmd)
  except OSError:
    raise Exception('Command "%s" does not exist' % cmd)
