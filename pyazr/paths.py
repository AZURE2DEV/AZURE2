"""The rule for a path that a project names, as the engine has it (experiment/ProjectPaths.h)."""


def is_absolute(path):
    """Whether ``path`` is absolute: a leading ``/`` or ``\\``, or a drive letter with its separator (``C:/...``).

    The twin of the engine's ``IsAbsolutePath``, the same on every system.  ``os.path.isabs`` is the system's rule:
    on Windows ``/tmp/x`` is not absolute for it (Python 3.13 on), and on Linux ``C:/x`` is not."""
    path = str(path)
    if path[:1] in ("/", "\\"):
        return True
    return len(path) > 2 and path[0].isascii() and path[0].isalpha() and path[1] == ":" and path[2] in ("/", "\\")
