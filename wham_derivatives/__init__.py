__all__ = []

from pathlib import Path
from importlib import import_module
from sys import modules

package = modules[__name__]
initfile = Path(__file__)
for entry in initfile.parent.iterdir():
    is_file = entry.is_file()
    is_pyfile = entry.name.endswith('.py')
    is_not_initpy = (entry != initfile)
    if is_file and is_pyfile and is_not_initpy:
        module_name = entry.name.removesuffix('.py')
        module_path = __name__ + '.' + module_name
        module = import_module(module_path)
        setattr(package, module_name, module)
        __all__.append(module_name)