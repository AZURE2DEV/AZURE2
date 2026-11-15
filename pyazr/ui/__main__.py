"""``python -m pyazr.ui project.azr``: serve a project on 127.0.0.1 (see :func:`pyazr.ui.server.main`)."""
import argparse


def cli(argv=None):
    ap = argparse.ArgumentParser(prog="python -m pyazr.ui", description="Serve an AZURE2 project to the web UI's API "
                                 "on a loopback address; prints the URL with its token.")
    ap.add_argument("project", help="the format-3 project (.azr); a new .azr in an existing folder starts an empty "
                    "project there, written when saved")
    ap.add_argument("--port", type=int, default=0, help="the TCP port (default: a free one above 20000)")
    ap.add_argument("--threads", type=int, default=1, help="the compute budget in threads (default 1)")
    ap.add_argument("--token", default=None, help="a fixed access token (default: a random one)")
    ap.add_argument("--host", default="127.0.0.1", help="a loopback address (anything else is refused)")
    ap.add_argument("--root", default=None, help="the directory files may be read and written in (default: the "
                    "project's)")
    ap.add_argument("--log-requests", action="store_true", help="log every request (method, path, status) on "
                    "standard error, as uvicorn's access log")
    a = ap.parse_args(argv)
    try:
        from .server import main
    except ImportError as e:
        raise SystemExit(f"pyazr.ui needs fastapi and uvicorn ({e}); pip install pyazr[web]")
    main(a.project, a.port, a.threads, a.token, host=a.host, root=a.root, log_requests=a.log_requests)


if __name__ == "__main__":
    cli()
