"""
Test for metareference APIs.
"""
from django.contrib.auth import get_user_model
from django.test import TestCase
from django.urls import reverse
from rest_framework import status
from rest_framework.test import APIClient
from core.models import (Linelist, Species, SpeciesMetadata,
                         Reference, MetaReference)
from data.serializers import MetaReferenceSerializer
import json
from rdkit import Chem
from rdkit.Chem import Descriptors
import selfies as sf


def create_linelist(linelist_name='Test Linelist'):
    """Helper function to create a linelist."""
    return Linelist.objects.create(linelist_name=linelist_name)


def create_species(**params):
    """Helper function to create a species."""
    defaults = {
        'name': json.dumps(['common_name', 'Test Species']),
        'iupac_name': 'Test IUPAC Name',
        'name_formula': 'Test Name Formula',
        'name_html': 'Test Name HTML',
        'molecular_mass': Descriptors.ExactMolWt(Chem.MolFromSmiles('CC')),
        'smiles': 'CC',
        'standard_inchi': 'test inchi',
        'standard_inchi_key': 'test inchi',
        'selfies': sf.encoder('CC'),
        'mol_obj': 'CC',
        'notes': 'Test Species',
    }
    defaults.update(params)

    return Species.objects.create(**defaults)


def create_meta(species_id, linelist_id, **params):
    """Helper function to create species metadata."""
    defaults = {
        'species_id': species_id,
        'molecule_tag': 1,
        'hyperfine': False,
        'degree_of_freedom': 3,
        'category': 'asymmetric top',
        'partition_function': json.dumps({'300.000': '331777.6674'}),
        'mu_a': 0.5,
        'mu_b': 0,
        'mu_c': 0,
        'a_const': 0.123,
        'b_const': 0,
        'c_const': 0,
        'linelist_id': linelist_id,
        'data_date': '2020-01-01',
        'data_contributor': 'Test Contributor',
        'int_file': 'test_int_file',
        'var_file': 'test_var_file',
        'fit_file': 'test_fit_file',
        'lin_file': 'test_lin_file',
        'qpart_file': 'test_qpart_file',
        'notes': 'Test Species Metadata',
    }
    defaults.update(params)

    return SpeciesMetadata.objects.create(**defaults)


def create_reference(**params):
    """Helper function to create a reference."""
    defaults = {
        'doi': '10.1021/acs.jcim.0c00128',
        'ref_url': 'https://doi.org/10.1021/acs.jcim.0c00128',
        'bibtex': 'bibtex_file',
        'notes': 'Test reference',
    }
    defaults.update(params)

    return Reference.objects.create(**defaults)


def create_metaref(meta_id, ref_id, **params):
    """Helper function to create a metareference."""
    defaults = {
        'meta_id': meta_id,
        'ref_id': ref_id,
        'dipole_moment': True,
        'spectrum': False,
        'notes': 'Test reference',
    }
    defaults.update(params)

    return MetaReference.objects.create(**defaults)


class PublicMetarefApiTests(TestCase):
    """Test the publicly available metareference API."""

    def setUp(self):
        self.client = APIClient()

    def test_auth_required_for_post(self):
        """Test that authentication is required for post
        creating metareference."""
        url = reverse('data:metareference-list')
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        ref = create_reference()
        payload = {
            'meta': meta.id,
            'ref': ref.id,
            'dipole_moment': True,
            'spectrum': False,
            'notes': 'Test reference'}

        res = self.client.post(url, payload)
        self.assertEqual(res.status_code, status.HTTP_401_UNAUTHORIZED)

    def test_auth_required_for_put(self):
        """Test that authentication is required for put updating
        metareference."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        ref = create_reference()
        metaref = create_metaref(meta.id, ref.id)
        url = reverse('data:metareference-detail', args=[metaref.id])
        payload = {
            'meta': meta.id,
            'ref': ref.id,
            'dipole_moment': True,
            'spectrum': True,
            'notes': 'Test put',
            '_change_reason': 'Test change reason'}
        res = self.client.put(url, payload)
        self.assertEqual(res.status_code, status.HTTP_401_UNAUTHORIZED)

    def test_auth_required_for_patch(self):
        """Test that authentication is required for patch updating
        metareference."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        ref = create_reference()
        metaref = create_metaref(meta.id, ref.id)
        url = reverse('data:metareference-detail', args=[metaref.id])
        payload = {'notes': 'test patch',
                   '_change_reason': 'Test change reason'}
        res = self.client.patch(url, payload)
        self.assertEqual(res.status_code, status.HTTP_401_UNAUTHORIZED)

    def test_auth_required_for_delete(self):
        """Test that authentication is required for deleting metareference."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        ref = create_reference()
        metaref = create_metaref(meta.id, ref.id)
        url = reverse('data:metareference-detail',
                      args=[metaref.id]) + '?delete_reason=Test delete reason'
        res = self.client.delete(url)
        self.assertEqual(res.status_code, status.HTTP_401_UNAUTHORIZED)

    def test_get_metaref_list(self):
        """Test getting a list of metareference."""
        species = create_species()
        linelist1 = create_linelist(linelist_name='Test Linelist 1')
        linelist2 = create_linelist(linelist_name='Test Linelist 2')
        meta1 = create_meta(species.id, linelist1.id)
        meta2 = create_meta(species.id, linelist2.id)
        ref = create_reference()
        create_metaref(meta1.id, ref.id)
        create_metaref(meta2.id, ref.id)
        url = reverse('data:metareference-list')

        res = self.client.get(url)
        metarefs = MetaReference.objects.all().order_by('-id')
        serializer = MetaReferenceSerializer(metarefs, many=True)
        self.assertEqual(res.status_code, status.HTTP_200_OK)
        self.assertEqual(res.data, serializer.data)

    def test_get_metaref_detail(self):
        """Test getting a metareference detail."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        ref = create_reference()
        metaref = create_metaref(meta.id, ref.id)
        url = reverse('data:metareference-detail', args=[metaref.id])

        res = self.client.get(url)
        serializer = MetaReferenceSerializer(metaref)
        self.assertEqual(res.status_code, status.HTTP_200_OK)
        self.assertEqual(res.data, serializer.data)


class PrivateMetarefApiTests(TestCase):
    """Test the private metareference API."""

    def setUp(self):
        self.client = APIClient()
        self.user = get_user_model().objects.create_user(
            email='example@example.com',
            password='testpass',
            name='Test User',
            organization='Test Organization')
        self.client.force_authenticate(self.user)

    def test_create_metaref(self):
        """Test creating a metaference."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        ref = create_reference()
        payload = {'meta': meta.id,
                   'ref': ref.id,
                   'dipole_moment': True,
                   'spectrum': False,
                   'notes': 'Test reference'}
        url = reverse('data:metareference-list')
        res = self.client.post(url, payload)
        self.assertEqual(res.status_code, status.HTTP_201_CREATED)
        metaref = MetaReference.objects.get(id=res.data['id'])
        history_exists = MetaReference.history.filter(id=metaref.id).exists()
        self.assertTrue(history_exists)

    def test_partial_update(self):
        """Test updating a metareference with patch."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        ref = create_reference()
        metaref = create_metaref(meta.id, ref.id)
        url = reverse('data:metareference-detail', args=[metaref.id])
        payload = {'notes': 'Updated metaref',
                   '_change_reason': 'Test patch'}
        res = self.client.patch(url, payload)
        self.assertEqual(res.status_code, status.HTTP_200_OK)
        linelist.refresh_from_db()
        self.assertEqual(MetaReference.history.filter(id=metaref.id).first(
        ).history_change_reason, payload['_change_reason'])
        self.assertEqual(MetaReference.history.filter(id=metaref.id).first(
        ).history_user_id, self.user.id)
        history_count = MetaReference.history.filter(id=metaref.id).count()
        self.assertEqual(history_count, 2)

    def test_full_update(self):
        """Test updating a metareference with put."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        ref = create_reference()
        metaref = create_metaref(meta.id, ref.id)
        url = reverse('data:metareference-detail', args=[metaref.id])
        payload = {'meta': meta.id,
                   'ref': ref.id,
                   'dipole_moment': True,
                   'spectrum': True,
                   'notes': 'Updated metaref put',
                   '_change_reason': 'Test put metaref'}
        res = self.client.put(url, payload)
        self.assertEqual(res.status_code, status.HTTP_200_OK)
        linelist.refresh_from_db()
        self.assertEqual(MetaReference.history.filter(id=metaref.id).first(
        ).history_change_reason, payload['_change_reason'])
        self.assertEqual(MetaReference.history.filter(id=metaref.id).first(
        ).history_user_id, self.user.id)
        history_count = MetaReference.history.filter(id=metaref.id).count()
        self.assertEqual(history_count, 2)

    def test_full_update_metaref_without_reason_fails(self):
        """Test updating a metareference with put without change
        reason fails."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        ref = create_reference()
        metaref = create_metaref(meta.id, ref.id)
        url = reverse('data:metareference-detail', args=[metaref.id])
        payload = {'meta': meta.id,
                   'ref': ref.id,
                   'dipole_moment': True,
                   'spectrum': True,
                   'notes': 'Updated metaref put'}
        res = self.client.put(url, payload)
        self.assertEqual(res.status_code, status.HTTP_400_BAD_REQUEST)

    def test_delete_metaref(self):
        """Test deleting a metareference."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        ref = create_reference()
        metaref = create_metaref(meta.id, ref.id)
        url = reverse('data:metareference-detail',
                      args=[metaref.id]) + \
            '?delete_reason=Test delete reason metaref'
        res = self.client.delete(url)
        self.assertEqual(res.status_code, status.HTTP_204_NO_CONTENT)
        self.assertFalse(MetaReference.objects.filter(id=metaref.id).exists())
        self.assertEqual(MetaReference.history.filter(id=metaref.id).first(
        ).history_change_reason, 'Test delete reason metaref')
        self.assertEqual(MetaReference.history.filter(id=metaref.id).first(
        ).history_user_id, self.user.id)
        history_count = MetaReference.history.filter(id=metaref.id).count()
        self.assertEqual(history_count, 2)

    def test_delete_metaref_without_delete_reason_fails(self):
        """Test deleting a metareference without delete reason fails."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        ref = create_reference()
        metaref = create_metaref(meta.id, ref.id)
        url = reverse('data:metareference-detail', args=[metaref.id])
        res = self.client.delete(url)
        self.assertEqual(res.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertTrue(MetaReference.objects.filter(id=metaref.id).exists())
