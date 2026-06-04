from .client import client, ClientError
from .server import server, ServerError
from .azure2 import azure2

__all__ = ["azure2", "client", "server", "ClientError", "ServerError"]
__version__ = '2.0.0'
