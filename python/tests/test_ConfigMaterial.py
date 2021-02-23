import os

import pytest
import numpy as np

from damask import ConfigMaterial
from damask import Table
from damask import Rotation

@pytest.fixture
def ref_path(ref_path_base):
    """Directory containing reference results."""
    return ref_path_base/'ConfigMaterial'


class TestConfigMaterial:

    @pytest.mark.parametrize('fname',[None,'test.yaml'])
    def test_load_save(self,ref_path,tmp_path,fname):
        reference = ConfigMaterial.load(ref_path/'material.yaml')
        os.chdir(tmp_path)
        if fname is None:
            reference.save()
            new = ConfigMaterial.load('material.yaml')
        else:
            reference.save(fname)
            new = ConfigMaterial.load(fname)
        assert reference == new

    def test_valid_complete(self,ref_path):
        material_config = ConfigMaterial.load(ref_path/'material.yaml')
        assert material_config.is_valid and material_config.is_complete

    def test_invalid_lattice(self,ref_path):
        material_config = ConfigMaterial.load(ref_path/'material.yaml')
        material_config['phase']['Aluminum']['lattice']='fxc'
        assert not material_config.is_valid

    def test_invalid_orientation(self,ref_path):
        material_config = ConfigMaterial.load(ref_path/'material.yaml')
        material_config['material'][0]['constituents'][0]['O']=[0,0,0,0]
        assert not material_config.is_valid

    def test_invalid_fraction(self,ref_path):
        material_config = ConfigMaterial.load(ref_path/'material.yaml')
        material_config['material'][0]['constituents'][0]['v']=.9
        assert not material_config.is_valid

    @pytest.mark.parametrize('item',['homogenization','phase','material'])
    def test_incomplete_missing(self,ref_path,item):
        material_config = ConfigMaterial.load(ref_path/'material.yaml')
        del material_config[item]
        assert not material_config.is_complete

    @pytest.mark.parametrize('item',['O','phase'])
    def test_incomplete_material_constituent(self,ref_path,item):
        material_config = ConfigMaterial.load(ref_path/'material.yaml')
        del material_config['material'][0]['constituents'][0][item]
        assert not material_config.is_complete

    def test_incomplete_material_homogenization(self,ref_path):
        material_config = ConfigMaterial.load(ref_path/'material.yaml')
        del material_config['material'][0]['homogenization']
        assert not material_config.is_complete

    def test_incomplete_homogenization_N_constituents(self,ref_path):
        material_config = ConfigMaterial.load(ref_path/'material.yaml')
        for h in material_config['homogenization'].keys():
            del material_config['homogenization'][h]['N_constituents']
        assert not material_config.is_complete

    def test_incomplete_phase_lattice(self,ref_path):
        material_config = ConfigMaterial.load(ref_path/'material.yaml')
        del material_config['phase']['Aluminum']['lattice']
        assert not material_config.is_complete

    def test_incomplete_wrong_phase(self,ref_path):
        material_config = ConfigMaterial.load(ref_path/'material.yaml')
        new = material_config.material_rename_phase({'Steel':'FeNbC'})
        assert not new.is_complete

    def test_incomplete_wrong_homogenization(self,ref_path):
        material_config = ConfigMaterial.load(ref_path/'material.yaml')
        new = material_config.material_rename_homogenization({'Taylor':'isostrain'})
        assert not new.is_complete

    def test_from_table(self):
        N = np.random.randint(3,10)
        a = np.vstack((np.hstack((np.arange(N),np.arange(N)[::-1])),np.ones(N*2),np.zeros(N*2),np.ones(N*2),np.ones(N*2))).T
        print(a)
        t = Table(a,{'varying':1,'constant':4})
        c = ConfigMaterial.from_table(t,**{'phase':'varying','O':'constant','homogenization':'4_constant'})
        assert len(c['material']) == N
        for i,m in enumerate(c['material']):
            assert m['homogenization'] == 1 and (m['constituents'][0]['O'] == [1,0,1,1]).all()

    @pytest.mark.parametrize('N,n,kw',[
                                        (1,1,{'phase':'Gold',
                                              'O':[1,0,0,0],
                                              'homogenization':'SX'}),
                                        (3,1,{'phase':'Gold',
                                              'O':Rotation.from_random(3),
                                              'homogenization':'SX'}),
                                        (2,3,{'phase':np.broadcast_to(['a','b','c'],(2,3)),
                                              'O':Rotation.from_random((2,3)),
                                              'homogenization':['SX','PX']}),
                                        ])
    def test_material_add(self,kw,N,n):
        m = ConfigMaterial().material_add(**kw)
        assert len(m['material']) == N
        assert len(m['material'][0]['constituents']) == n
