"""Write a file whole or not at all.

The GUI reads result files while analysis runs write them, so writing in
place lets a reader see half a file, and a write that fails part-way destroys
the last good one. ``atomic_open`` writes a temporary file in the same folder
and renames it over the target only once it is complete.

The temporary file is created like ``open()`` creates a file, so the user's
umask applies (``tempfile.mkstemp`` would make it 0600 and take group read
access away from shared output folders). Its name, ``<target>.<random>.tmp``,
matches no reader's exact name or ``*.npy`` / ``*.json`` glob.
"""

import contextlib
import os
import secrets
from typing import IO, Iterator


@contextlib.contextmanager
def atomic_open(path: str, mode: str = "w") -> Iterator[IO]:
    """Yield a file to write; on success it replaces ``path``, on error it is removed.

    ``mode`` is ``"w"`` (text) or ``"wb"`` (binary).
    """
    path = os.fspath(path)
    folder = os.path.dirname(os.path.abspath(path))
    base = os.path.basename(path)
    while True:
        tmp_path = os.path.join(folder, f"{base}.{secrets.token_hex(6)}.tmp")
        try:
            # 0o666 & ~umask, exactly as open() does
            fd = os.open(tmp_path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o666)
            break
        except FileExistsError:
            continue
    try:
        with os.fdopen(fd, mode) as f:
            yield f
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp_path, path)
    except BaseException:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass
        raise
