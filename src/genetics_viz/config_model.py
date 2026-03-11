"""Application configuration model — YAML config loading, saving, and user management."""

import fcntl
import hashlib
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
    password: str  # sha512 hex digest
    role: str  # "reader" | "curator" | "administrator"


@dataclass
class AppConfig:
    """Top-level application configuration."""

    config_path: Path
    data_directories: list[DataDirectoryConfig] = field(default_factory=list)
    user_list: list[UserConfig] = field(default_factory=list)
    storage_secret: str = ""


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


def hash_password(plain: str) -> str:
    """Hash a plaintext password to a sha512 hex digest."""
    return hashlib.sha512(plain.encode()).hexdigest()


def verify_password(plain: str, hashed: str) -> bool:
    """Verify a plaintext password against a sha512 hex digest."""
    return hashlib.sha512(plain.encode()).hexdigest() == hashed


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
    for entry in raw.get("data_directories", []):
        dirs.append(
            DataDirectoryConfig(
                path=str(entry.get("path", "")),
                description=str(entry.get("description", "")),
                default=bool(entry.get("default", False)),
            )
        )

    # Parse users
    users: list[UserConfig] = []
    for entry in raw.get("user_list", []):
        users.append(
            UserConfig(
                username=str(entry.get("username", "")),
                password=str(entry.get("password", "")),
                role=str(entry.get("role", "reader")),
            )
        )

    # Storage secret — auto-generate if missing
    storage_secret = raw.get("storage_secret", "")
    secret_was_missing = not storage_secret
    if secret_was_missing:
        storage_secret = secrets.token_hex(32)

    config = AppConfig(
        config_path=config_path,
        data_directories=dirs,
        user_list=users,
        storage_secret=storage_secret,
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
