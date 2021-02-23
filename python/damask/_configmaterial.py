import numpy as np

from . import Config
from . import Rotation
from . import Orientation
from . import util

class ConfigMaterial(Config):
    """Material configuration."""

    _defaults = {'material': [],
                 'homogenization': {},
                 'phase': {}}

    def __init__(self,d=_defaults):
        """Initialize object with default dictionary keys."""
        super().__init__(d)


    def save(self,fname='material.yaml',**kwargs):
        """
        Save to yaml file.

        Parameters
        ----------
        fname : file, str, or pathlib.Path, optional
            Filename or file for writing. Defaults to 'material.yaml'.
        **kwargs
            Keyword arguments parsed to yaml.dump.

        """
        super().save(fname,**kwargs)


    @classmethod
    def load(cls,fname='material.yaml'):
        """
        Load from yaml file.

        Parameters
        ----------
        fname : file, str, or pathlib.Path, optional
            Filename or file for writing. Defaults to 'material.yaml'.

        """
        return super(ConfigMaterial,cls).load(fname)


    @staticmethod
    def from_table(table,**kwargs):
        """
        Load from an ASCII table.

        Parameters
        ----------
        table : damask.Table
            Table that contains material information.
        **kwargs
            Keyword arguments where the key is the name and the value specifies
            the label of the data column in the table.

        Examples
        --------
        >>> import damask
        >>> import damask.ConfigMaterial as cm
        >>> t = damask.Table.load('small.txt')
        >>> t
            pos  pos  pos   qu   qu    qu    qu   phase    homog
        0    0    0    0  0.19  0.8   0.24 -0.51  Aluminum SX
        1    1    0    0  0.8   0.19  0.24 -0.51  Steel    SX
        1    1    1    0  0.8   0.19  0.24 -0.51  Steel    SX
        >>> cm.from_table(t,O='qu',phase='phase',homogenization='homog')
        material:
          - constituents:
              - O: [0.19, 0.8, 0.24, -0.51]
                v: 1.0
                phase: Aluminum
            homogenization: SX
          - constituents:
              - O: [0.8, 0.19, 0.24, -0.51]
                v: 1.0
                phase: Steel
            homogenization: SX
        homogenization: {}
        phase: {}

        """
        kwargs_       = {k:table.get(v) for k,v in kwargs.items()}

        _,idx = np.unique(np.hstack(list(kwargs_.values())),return_index=True,axis=0)
        idx = np.sort(idx)
        kwargs_ = {k:np.atleast_1d(v[idx].squeeze()) for k,v in kwargs_.items()}

        return ConfigMaterial().material_add(**kwargs_)


    @property
    def is_complete(self):
        """Check for completeness."""
        ok = True
        for top_level in ['homogenization','phase','material']:
            ok &= top_level in self
            if top_level not in self: print(f'{top_level} entry missing')

        if ok:
           ok &= len(self['material']) > 0
           if len(self['material']) < 1: print('Incomplete material definition')

        if ok:
            homogenization = set()
            phase          = set()
            for i,v in enumerate(self['material']):
                if 'homogenization' in v:
                    homogenization.add(v['homogenization'])
                else:
                    print(f'No homogenization specified in material {i}')
                    ok = False

                if 'constituents' in v:
                    for ii,vv in enumerate(v['constituents']):
                        if 'O' not in vv:
                            print('No orientation specified in constituent {ii} of material {i}')
                            ok = False
                        if 'phase' in vv:
                            phase.add(vv['phase'])
                        else:
                            print(f'No phase specified in constituent {ii} of material {i}')
                            ok = False

            for k,v in self['phase'].items():
                if 'lattice' not in v:
                    print(f'No lattice specified in phase {k}')
                    ok = False

            for k,v in self['homogenization'].items():
                if 'N_constituents' not in v:
                    print(f'No. of constituents not specified in homogenization {k}')
                    ok = False

            if phase - set(self['phase']):
                print(f'Phase(s) {phase-set(self["phase"])} missing')
                ok = False
            if homogenization - set(self['homogenization']):
                print(f'Homogenization(s) {homogenization-set(self["homogenization"])} missing')
                ok = False

        return ok


    @property
    def is_valid(self):
        """Check for valid file layout."""
        ok = True

        if 'phase' in self:
            for k,v in self['phase'].items():
                if 'lattice' in v:
                    try:
                        Orientation(lattice=v['lattice'])
                    except KeyError:
                        s = v['lattice']
                        print(f"Invalid lattice: '{s}' in phase '{k}'")
                        ok = False

        if 'material' in self:
            for i,m in enumerate(self['material']):
                if 'constituents' in m:
                    v = 0.0
                    for c in m['constituents']:
                        v+= float(c['v'])
                        if 'O' in c:
                            try:
                                Rotation.from_quaternion(c['O'])
                            except ValueError:
                                o = c['O']
                                print(f"Invalid orientation: '{o}' in material '{i}'")
                                ok = False
                    if not np.isclose(v,1.0):
                        print(f"Invalid total fraction (v) '{v}' in material '{i}'")
                        ok = False

        return ok


    def material_rename_phase(self,mapping,ID=None,constituent=None):
        """
        Change phase name in material.

        Parameters
        ----------
        mapping: dictionary
            Mapping from old name to new name
        ID: list of ints, optional
            Limit renaming to selected material IDs.
        constituent: list of ints, optional
            Limit renaming to selected constituents.

        """
        dup = self.copy()
        for i,m in enumerate(dup['material']):
            if ID is not None and i not in ID: continue
            for c in m['constituents']:
                if constituent is not None and c not in constituent: continue
                try:
                    c['phase'] = mapping[c['phase']]
                except KeyError:
                    continue
        return dup


    def material_rename_homogenization(self,mapping,ID=None):
        """
        Change homogenization name in material.

        Parameters
        ----------
        mapping: dictionary
            Mapping from old name to new name
        ID: list of ints, optional
            Limit renaming to selected homogenization IDs.

        """
        dup = self.copy()
        for i,m in enumerate(dup['material']):
            if ID is not None and i not in ID: continue
            try:
                m['homogenization'] = mapping[m['homogenization']]
            except KeyError:
                continue
        return dup


    def material_add(self,**kwargs):
        """
        Add material entries.

        Parameters
        ----------
        constituents : dict, optional
            Entries for 'constituents' as key-value pair.
        **kwargs
            Key-value pairs.

        Examples
        --------
        >>> import damask
        >>> m = damask.ConfigMaterial().material_add(phase          = ['Aluminum','Steel','Aluminum'],
        ...                                          O              = damask.Rotation.from_random(3)},
        ...                                          homogenization = 'SX')
        >>> m
        material:
          - constituents:
              - O: [0.577764, -0.146299, -0.617669, 0.513010]
                v: 1.0
                phase: Aluminum
            homogenization: SX
          - constituents:
              - O: [0.184176, 0.340305, 0.737247, 0.553840]
                v: 1.0
                phase: Steel
            homogenization: SX
          - constituents:
              - O: [0.0886257, -0.144848, 0.615674, -0.769487]
                v: 1.0
                phase: Aluminum
            homogenization: SX
        homogenization: {}
        phase: {}

        """
        N,n,shaped = 1,1,{}

        for k,v in kwargs.items():
            shaped[k] = np.array(v)
            s = shaped[k].shape[:-1] if k=='O' else shaped[k].shape
            N = max(N,s[0]) if len(s)>0 else N
            n = max(n,s[1]) if len(s)>1 else n

        mat = [{'constituents':[{} for _ in range(n)]} for _ in range(N)]

        if 'v' not in kwargs:
            shaped['v'] = np.broadcast_to(1/n,(N,n))

        for k,v in shaped.items():
            obj = np.broadcast_to(v.reshape(util.shapeshifter(v.shape,(N,n,4))),(N,n,4)) if k=='O' else \
                  np.broadcast_to(v.reshape(util.shapeshifter(v.shape,(N,n  ))),(N,n))
            for i in range(N):
                if k in ['phase','O','v']:
                    for j in range(n):
                        mat[i]['constituents'][j][k] = obj[i,j].item() if isinstance(obj[i,j],np.generic) else obj[i,j]
                else:
                    mat[i][k] = obj[i,0].item() if isinstance(obj[i,0],np.generic) else obj[i,0]

        dup = self.copy()
        dup['material'] = dup['material'] + mat if 'material' in dup else mat

        return dup
