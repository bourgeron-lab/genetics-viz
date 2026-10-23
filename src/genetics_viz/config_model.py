"""Application configuration model — YAML config loading, saving, and user management."""

import fcntl
import hashlib
import hmac
import secrets
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import yaml

# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------


@dataclass
class DataDirectoryConfig:
    """A single data directory entry from the YAML config."""

    path: str
    description: str = ""
    default: bool = False


@dataclass
class UserConfig:
    """A single user entry from the YAML config."""

    username: str
    password: str  # scrypt hash (legacy: unsalted sha512 hex digest)
    role: str  # "reader" | "curator" | "administrator"


@dataclass
class AppConfig:
    """Top-level application configuration."""

    config_path: Path
    data_directories: list[DataDirectoryConfig] = field(default_factory=list)
    user_list: list[UserConfig] = field(default_factory=list)
    storage_secret: str = ""
    poll_interval: int = 30
    ancestry_reference_dir: str = ""
    reference_fasta: str = ""
    genetrek_url: str = ""
    genetrek_cache_dir: str = ""
    gene_index_dir: str = ""


#: GeneTrek instance used when the config names none. The dev instance, until
#: the production one serves the same API.
DEFAULT_GENETREK_URL = "https://genetrek-py.dev.pasteur.cloud/"


# ---------------------------------------------------------------------------
# Module-level singleton
# ---------------------------------------------------------------------------

_app_config: Optional[AppConfig] = None


def get_config() -> AppConfig:
    """Return the global AppConfig. Raises RuntimeError if not loaded."""
    if _app_config is None:
        raise RuntimeError("AppConfig not loaded — call load_config() first")
    return _app_config


def get_config_path() -> Path:
    """Return the path to the YAML config file."""
    return get_config().config_path


# ---------------------------------------------------------------------------
# Password helpers
# ---------------------------------------------------------------------------


#: scrypt cost: N = 2**15, r = 8 takes 32 MiB and about 50 ms per hash. The
#: parameters are stored in each hash, so raising them later leaves existing
#: hashes verifiable (and ``needs_rehash`` upgrades them at the next login).
SCRYPT_N = 2**15
SCRYPT_R = 8
SCRYPT_P = 1
_SCRYPT_MAXMEM = 64 * 2**20
_SCRYPT_PREFIX = "scrypt$"


def _scrypt(plain: str, salt: bytes, n: int, r: int, p: int) -> bytes:
    return hashlib.scrypt(
        plain.encode(), salt=salt, n=n, r=r, p=p, maxmem=_SCRYPT_MAXMEM, dklen=64
    )


def hash_password(plain: str) -> str:
    """Hash a plaintext password with salted scrypt.

    Format: ``scrypt$<N>$<r>$<p>$<salt hex>$<hash hex>``.
    """
    salt = secrets.token_bytes(16)
    digest = _scrypt(plain, salt, SCRYPT_N, SCRYPT_R, SCRYPT_P)
    return (
        f"{_SCRYPT_PREFIX}{SCRYPT_N}${SCRYPT_R}${SCRYPT_P}${salt.hex()}${digest.hex()}"
    )


def verify_password(plain: str, hashed: str) -> bool:
    """Verify a plaintext password against a scrypt hash or a legacy SHA-512 one.

    Legacy hashes are unsalted SHA-512 hex digests, from before scrypt. They
    still verify, so nobody needs a reset; ``needs_rehash`` tells the login to
    replace them.
    """
    if hashed.startswith(_SCRYPT_PREFIX):
        try:
            n, r, p, salt, digest = hashed[len(_SCRYPT_PREFIX) :].split("$")
            computed = _scrypt(plain, bytes.fromhex(salt), int(n), int(r), int(p))
        except ValueError:
            return False
        return hmac.compare_digest(computed.hex(), digest)
    legacy = hashlib.sha512(plain.encode()).hexdigest()
    return bool(hashed) and hmac.compare_digest(legacy, hashed)


def needs_rehash(hashed: str) -> bool:
    """Return True for a legacy SHA-512 hash or scrypt with older parameters."""
    current = f"{_SCRYPT_PREFIX}{SCRYPT_N}${SCRYPT_R}${SCRYPT_P}$"
    return not hashed.startswith(current)


def generate_random_password() -> str:
    """Generate a random URL-safe password (16 chars)."""
    return secrets.token_urlsafe(12)


# ---------------------------------------------------------------------------
# Load / save
# ---------------------------------------------------------------------------


def load_config(config_path: Path) -> AppConfig:
    """Load the YAML configuration file and set the global singleton.

    If ``storage_secret`` is missing from the file, a random one is generated
    and written back so that sessions persist across restarts.
    """
    global _app_config

    with open(config_path, "r") as f:
        raw = yaml.safe_load(f) or {}

    # Parse data directories
    dirs: list[DataDirectoryConfig] = []
    # ``or`` rather than a .get default: a key present but empty loads as None.
    for entry in raw.get("data_directories") or []:
        dirs.append(
            DataDirectoryConfig(
                path=str(entry.get("path") or ""),
                description=str(entry.get("description") or ""),
                default=bool(entry.get("default")),
            )
        )

    # Parse users
    users: list[UserConfig] = []
    for entry in raw.get("user_list") or []:
        users.append(
            UserConfig(
                username=str(entry.get("username") or ""),
                password=str(entry.get("password") or ""),
                role=str(entry.get("role") or "reader"),
            )
        )

    # Storage secret — auto-generate if missing
    storage_secret = raw.get("storage_secret") or ""
    secret_was_missing = not storage_secret
    if secret_was_missing:
        storage_secret = secrets.token_hex(32)

    poll_interval = int(raw.get("poll_interval") or 30)
    ancestry_reference_dir = str(raw.get("ancestry_reference_dir") or "")
    reference_fasta = str(raw.get("reference_fasta") or "")
    genetrek_url = str(raw.get("genetrek_url") or "")
    genetrek_cache_dir = str(raw.get("genetrek_cache_dir") or "")
    gene_index_dir = str(raw.get("gene_index_dir") or "")

    config = AppConfig(
        config_path=config_path,
        data_directories=dirs,
        user_list=users,
        storage_secret=storage_secret,
        poll_interval=poll_interval,
        ancestry_reference_dir=ancestry_reference_dir,
        reference_fasta=reference_fasta,
        genetrek_url=genetrek_url,
        genetrek_cache_dir=genetrek_cache_dir,
        gene_index_dir=gene_index_dir,
    )
    _app_config = config

    # Persist the generated secret so it survives restarts
    if secret_was_missing:
        save_config(config)

    return config


def save_config(config: AppConfig) -> None:
    """Write the current config back to the YAML file with file locking."""
    data: dict = {
        "data_directories": [
            {
                k: v
                for k, v in {
                    "path": d.path,
                    "description": d.description,
                    "default": d.default if d.default else None,
                }.items()
                if v is not None and v != ""
            }
            for d in config.data_directories
        ],
        "user_list": [
            {
                "username": u.username,
                "password": u.password,
                "role": u.role,
            }
            for u in config.user_list
        ],
        "storage_secret": config.storage_secret,
    }

    # Emit the remaining top-level keys only when they carry information.
    # ``save_config`` rewrites the whole file from known fields, so anything
    # omitted here is silently dropped the next time an admin saves.
    if config.ancestry_reference_dir:
        data["ancestry_reference_dir"] = config.ancestry_reference_dir
    if config.reference_fasta:
        data["reference_fasta"] = config.reference_fasta
    if config.poll_interval != 30:
        data["poll_interval"] = config.poll_interval
    if config.genetrek_url:
        data["genetrek_url"] = config.genetrek_url
    if config.genetrek_cache_dir:
        data["genetrek_cache_dir"] = config.genetrek_cache_dir
    if config.gene_index_dir:
        data["gene_index_dir"] = config.gene_index_dir

    with open(config.config_path, "w") as f:
        fcntl.flock(f.fileno(), fcntl.LOCK_EX)
        try:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)
        finally:
            fcntl.flock(f.fileno(), fcntl.LOCK_UN)


def get_default_data_dir(config: AppConfig) -> str:
    """Return the path string of the default data directory."""
    for d in config.data_directories:
        if d.default:
            return d.path
    # Fallback: first directory
    if config.data_directories:
        return config.data_directories[0].path
    return ""


def get_ancestry_reference_dir() -> Optional[Path]:
    """Return the ancestry/PGS reference bundle root, or ``None`` when unset.

    The directory holds one ``v<bundle_version>/`` subdirectory per reference
    bundle; the version to use is derived from the per-family data filenames
    rather than configured here.
    """
    if _app_config is None:
        return None
    if not _app_config.ancestry_reference_dir:
        return None
    return Path(_app_config.ancestry_reference_dir)


def get_reference_fasta() -> Optional[Path]:
    """Return the FASTA the CRAMs were aligned to, or ``None`` when unset.

    pysam needs it to decode CRAM reads for the SV read-evidence panel. It must
    be the alignment reference itself (GRCh38_GIABv3), not a generic hg38.
    """
    if _app_config is None or not _app_config.reference_fasta:
        return None
    return Path(_app_config.reference_fasta)


def get_genetrek_url() -> str:
    """Return the GeneTrek base URL, always ending in ``/``."""
    url = (_app_config.genetrek_url if _app_config else "") or DEFAULT_GENETREK_URL
    return url if url.endswith("/") else url + "/"


def get_genetrek_cache_dir() -> Optional[Path]:
    """Return where the GeneTrek snapshot is cached, or ``None`` with no config.

    Defaults to ``.genetrek/`` next to the config file: that directory is
    already writable, since ``load_config`` may write the file back.
    """
    if _app_config is None:
        return None
    if _app_config.genetrek_cache_dir:
        return Path(_app_config.genetrek_cache_dir)
    return _app_config.config_path.parent / ".genetrek"


def get_gene_index_dir() -> Optional[Path]:
    """Return where the cross-cohort gene index is kept, or ``None`` with no config.

    Defaults to ``.gene_index/`` next to the config file, as for the GeneTrek
    cache: local disk, so lookups never cross the data mount. It holds
    genotypes, so it is created readable by the app's user only.
    """
    if _app_config is None:
        return None
    if _app_config.gene_index_dir:
        return Path(_app_config.gene_index_dir)
    return _app_config.config_path.parent / ".gene_index"


def get_preferences_path() -> Optional[Path]:
    """Return the per-user preferences file, next to the config file.

    Kept out of the config file itself: ``save_config`` rewrites that file
    whole from the in-memory singleton, so frequent profile saves there would
    race with admin saves.
    """
    if _app_config is None:
        return None
    path = _app_config.config_path
    return path.with_name(f"{path.stem}.preferences.yaml")
