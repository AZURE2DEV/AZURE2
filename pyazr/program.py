"""The program AZURE2 of an installed wheel.

A wheel carries the program inside the package, ``pyazr/bin/AZURE2``, and the console script ``AZURE2`` that pip
puts on the PATH is :func:`main`, which replaces itself by that program with the same arguments.  The program is not
in the wheel's scripts because a script is moved at installation: the libraries the wheel carries (GSL, the OpenMP
run time; copied in by auditwheel or delocate) are found by a path relative to the file that needs them, and that
path holds only for a file that stays inside the package.  In a checkout built with CMake there is no
``pyazr/bin``; the program is ``build/src/AZURE2`` and this module only says so.

On Windows the wheel is built with MinGW (the toolchain of the Windows program) and repaired by delvewheel, which
puts the DLLs the program loads (GSL, the OpenMP run time, the MinGW run time) in ``pyazr.libs`` beside the package,
unmangled; the program is ``pyazr/bin/AZURE2.exe``, :func:`main` runs it as a child process (Windows has no exec)
with that directory on its ``PATH``.  Not yet run on a Windows machine (``.github/workflows/wheels.yml`` is the first).
"""

import os
import subprocess
import sys


def path():
    """The program's file in this package, or ``None`` when the package has none (a checkout, not a wheel)."""
    name = "AZURE2.exe" if os.name == "nt" else "AZURE2"
    file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "bin", name)
    return file if os.path.isfile(file) else None


def _libs_directory():
    """The directory delvewheel puts a Windows wheel's DLLs in (``pyazr.libs``, beside the package), or ``None``."""
    folder = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "pyazr.libs")
    return folder if os.path.isdir(folder) else None


def main(argv=None):
    """Run the package's program with the command line's arguments, in this process's place."""
    program = path()
    if program is None:
        sys.stderr.write("AZURE2: this pyazr holds no program (pyazr/bin/AZURE2 is in a wheel only); in a checkout "
                         "run the CMake build's src/AZURE2.\n")
        return 2
    arguments = [program] + list(sys.argv[1:] if argv is None else argv)
    sys.stdout.flush()
    sys.stderr.flush()
    if os.name == "nt":
        environment = dict(os.environ)
        libs = _libs_directory()
        if libs:
            environment["PATH"] = libs + os.pathsep + environment.get("PATH", "")
        return subprocess.call(arguments, env=environment)
    os.execv(program, arguments)
