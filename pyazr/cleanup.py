"""Find and reap AZURE2 instances nothing owns any more.

A pyazr session stops the processes it started, and since it registers an
:mod:`atexit` hook it does so even when a script forgets to call ``close()``.
What it cannot cover is the interpreter dying without running its handlers -- a
notebook kernel restarted, a ``kill -9``, a crash. Those leave AZURE2 processes
running with nobody to talk to them.

The distinction this module rests on is **parentage**, not name. A process
belonging to a script that is still running has that script as its parent. When
the parent dies the child is reparented to init (PID 1), and that -- not the
command line -- is what marks it as abandoned. So a concurrently running
script's instances are never candidates, however many of them there are.

From the command line::

    python -m pyazr.cleanup            # list what would be reaped
    python -m pyazr.cleanup --kill     # reap it
"""

import os
import signal
import subprocess
import sys
import time

__all__ = ["find_orphans", "kill_orphans"]


def _running_processes():
    """(pid, ppid, command) for the current user's processes."""
    if os.name != "posix":
        return []
    try:
        out = subprocess.run(
            ["ps", "-o", "pid=,ppid=,command=", "-u", str(os.getuid())],
            capture_output=True, text=True, timeout=15,
        ).stdout
    except (OSError, subprocess.SubprocessError):
        return []

    rows = []
    for line in out.splitlines():
        parts = line.split(None, 2)
        if len(parts) < 3:
            continue
        try:
            rows.append((int(parts[0]), int(parts[1]), parts[2]))
        except ValueError:
            continue
    return rows


def find_orphans():
    """AZURE2 API instances whose parent process is gone.

    Returns a list of ``(pid, command)``. Only processes started with
    ``--use-api`` are considered -- a GUI AZURE2 or an interactive console run
    is somebody's actual session, not litter -- and only those reparented to
    PID 1, which is what makes them unreachable.
    """
    orphans = []
    for pid, ppid, command in _running_processes():
        if "AZURE2" not in command or "--use-api" not in command:
            continue
        if ppid != 1:            # still owned by a live process; leave it alone
            continue
        if pid == os.getpid():
            continue
        orphans.append((pid, command))
    return orphans


def kill_orphans(dry_run=False, timeout=5.0):
    """Terminate abandoned AZURE2 API instances.

    SIGTERM first, then SIGKILL for anything still alive after ``timeout``.
    Returns the list of ``(pid, command)`` that were targeted; with
    ``dry_run=True`` nothing is signalled and the list is what *would* be.
    """
    orphans = find_orphans()
    if dry_run or not orphans:
        return orphans

    for pid, _ in orphans:
        try:
            os.kill(pid, signal.SIGTERM)
        except (ProcessLookupError, PermissionError, OSError):
            pass

    deadline = time.time() + timeout
    while time.time() < deadline:
        if not any(_alive(pid) for pid, _ in orphans):
            break
        time.sleep(0.1)

    for pid, _ in orphans:
        if _alive(pid):
            try:
                os.kill(pid, signal.SIGKILL)
            except (ProcessLookupError, PermissionError, OSError):
                pass
    return orphans


def _alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OSError:
        return False
    return True


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    do_kill = "--kill" in argv or "-k" in argv

    if os.name != "posix":
        print("pyazr.cleanup only supports POSIX systems.")
        return 1

    targets = kill_orphans(dry_run=not do_kill)
    if not targets:
        print("No abandoned AZURE2 API instances found.")
        return 0

    verb = "Reaped" if do_kill else "Would reap"
    print(f"{verb} {len(targets)} abandoned AZURE2 API instance(s):")
    for pid, command in targets:
        print(f"  {pid}  {command[:110]}")
    if not do_kill:
        print("\nRe-run with --kill to terminate them.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
