"""Which BAM is this? A BAM is known by its bytes, not by its path.

Each step that adds a BAM's data to a sample's totals keeps the keys of the BAMs it
has added, so that running it again for the same BAM (after a restart, or because
the start-up scan and the watcher both saw it, or for a copy of it in another
folder) adds nothing (JOB-2).
"""

import hashlib
import os

_CHUNK = 64 * 1024
_TAG = b"robin-bam-key-v1\0"


def bam_key(path: str) -> str:
    """The BAM's key: its size, its first 64 KiB and its last 64 KiB, hashed (blake2b-128).

    A BAM of 128 KiB or less is hashed whole. The same file at another path, or a copy
    of it, gives the same key. Two different BAMs are taken to differ in size or in
    those bytes: the last 64 KiB holds the last reads and the end-of-file block. This
    assumes a BAM is not changed once it is complete.
    """
    with open(path, "rb") as f:
        size = os.fstat(f.fileno()).st_size
        h = hashlib.blake2b(digest_size=16)
        h.update(_TAG)
        h.update(size.to_bytes(8, "little"))
        h.update(f.read(_CHUNK))
        if size > _CHUNK:
            f.seek(max(_CHUNK, size - _CHUNK))
            h.update(f.read(_CHUNK))
    return h.hexdigest()
