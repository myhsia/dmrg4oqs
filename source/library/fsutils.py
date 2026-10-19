"""Filesystem and IO helpers: config loading, output writers and loaders."""
import copy, os, yaml
import numpy as np
from pathlib import Path


def _apply_relative(value, base):
    if not isinstance(value, str) or base is None:
        return value
    s = value.strip()
    if len(s) < 2 or s[0] not in "*/+-":
        return value
    try:
        factor = float(s[1:])
    except ValueError:
        return value
    if s[0] == "*":
        return base * factor
    if s[0] == "/":
        return base / factor
    if s[0] == "+":
        return base + factor
    return base - factor


def _merge(base, override):
    merged = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _merge(merged[key], value)
        else:
            merged[key] = _apply_relative(value, merged.get(key))
    return merged


def load_config(prefix: str, config_dir: str = "./config"):
    with open(os.path.join(config_dir, prefix + ".yaml"), "r") as f:
        override = yaml.safe_load(f)
    base_name = override.pop("base", None)
    if base_name is None:
        return override
    base_name = os.path.splitext(str(base_name))[0]
    with open(os.path.join(config_dir, base_name + ".yaml"), "r") as f:
        params = yaml.safe_load(f)
    return _merge(params, override)


def ensure_dir(directory) -> None:
    """Create ``directory`` (and any missing parents) if it does not exist."""
    Path(directory).mkdir(parents = True, exist_ok = True)


def save_time(t, Y, directory, filename, identifier = "",
              compress: bool = False) -> None:
    filepath = Path(directory) / filename
    ensure_dir(filepath.parent)
    if compress:
        np.savez(filepath, t = np.asarray(t), y = np.asarray(Y),
                 identifier = np.array(identifier))
        return
    if np.iscomplexobj(Y):
        datout = np.column_stack((t, Y.real, Y.imag))
        header_list = (["t", "Re", "Im"] if identifier == "" else
                       ["t", f"Re[{identifier}]", f"Im[{identifier}]"])
    else:
        Y = np.asarray(Y)
        Y = Y[:, None] if Y.ndim == 1 else Y
        datout = np.column_stack((t, Y))
        prefix = identifier if identifier else "Y"
        header_list = ["t"] + [f"{prefix}{i}" for i in range(Y.shape[1])]
    fmt = "%15.6e" * datout.shape[1]
    header = f"{header_list[0]:>13}" + "".join(f"{h:>15}" for h in header_list[1:])
    np.savetxt(filepath, datout, fmt = fmt, header = header)


def load_npz(path):
    data = np.load(path)
    return data["t"], data["y"]


def moments_filename(case: str, projection: str = "matrix") -> str:
    """Canonical moments data filename ``moments_<case>.{mat,scl}.npz``,
    where ``projection`` is ``"matrix"`` (``mat``) or ``"scalar"`` (``scl``).
    """
    return f"moments_{case}.{'scl' if projection == 'scalar' else 'mat'}.npz"


def save_moments(path, arr) -> None:
    """Save a moments array under the ``moments`` key of an ``.npz`` archive.

    ``arr`` may be a scalar ``(n_moments,)`` or matrix ``(n_moments, d, d)``
    complex array; the native shape is preserved (no flattening).
    """
    filepath = Path(path)
    ensure_dir(filepath.parent)
    np.savez(filepath, moments = np.asarray(arr))


def load_moments(path) -> np.ndarray:
    """Load a moments array from an ``.npz`` archive (key ``moments``).

    Returns the native shape: ``(n_moments, d, d)`` for the matrix projection
    or ``(n_moments,)`` for the scalar projection.
    """
    with np.load(path) as f:
        return f["moments"]


def tab_str(s, N = 80, char = '=', border = '!'):
    """Pad ``s`` to a border line of constant width ``N`` for console logging."""
    n = len(s)
    n_left = (N - n - 4) // 2
    n_right = N - n - 4 - n_left
    return border + char * n_left + ' ' + s + ' ' + char * n_right + border
