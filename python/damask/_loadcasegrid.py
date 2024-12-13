from typing import Optional, Union, Dict, Any, List

from numpy import ma

from ._typehints import FileHandle
from ._yaml import NiceDumper
from . import YAML


class MaskedMatrixDumper(NiceDumper):
    """Format masked matrices."""

    def represent_data(self, data: Any):
        return super().represent_data(data.astype(object).filled('x') if isinstance(data, ma.core.MaskedArray) else data) # type: ignore[attr-defined]


class LoadcaseGrid(YAML):
    """Load case for grid solver."""

    def __init__(self,
                 config: Optional[Union[str,Dict[str,Any]]] = None,
                 *,
                 solver: Optional[Dict[str,str]] = None,
                 loadstep: Optional[List[Dict[str,Any]]] = None):
        """
        New grid solver load case.

        Parameters
        ----------
        config : dict or str, optional
            Grid solver load case. String needs to be valid YAML.
        solver : dict, optional
            Solver configuration.
            Defaults to an empty dict if 'config' is not given.
        loadstep : list of dict, optional
            Load step configuration.
            Defaults to an empty list if 'config' is not given.

        """
        kwargs: Dict[str,Union[Dict[str,str],List[Dict[str,Any]]]] = {}
        default: Union[List,Dict]
        for arg,value,default in [('solver',solver,{}),('loadstep',loadstep,[])]:                   # type: ignore[assignment]
            if value is not None:
                kwargs[arg] = value
            elif config is None:
                kwargs[arg] = default

        super().__init__(config,**kwargs)


    def save(self,
             fname: FileHandle,
             **kwargs):
        """
        Save to YAML file.

        Parameters
        ----------
        fname : file, str, or pathlib.Path
            Filename or file to write.
        **kwargs : dict
            Keyword arguments parsed to yaml.dump.

        """
        if 'Dumper' not in kwargs:
            kwargs['Dumper'] = MaskedMatrixDumper

        super().save(fname=fname,**kwargs)
