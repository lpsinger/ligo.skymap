import os
import pkgutil

__all__ = ()

# Import all symbols from all submodules of this module.
for _, module, _ in pkgutil.iter_modules([os.path.dirname(__file__)]):
    if module not in {"tests"}:
        exec(  # noqa: S102
            f"from . import {module};"
            f'__all__ += getattr({module}, "__all__", ());'
            f"from .{module} import *"
        )
    del module

# Clean up
del os, pkgutil
