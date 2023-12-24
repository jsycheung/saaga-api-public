"""
Test for reference APIs.
"""
from django.contrib.auth import get_user_model
from django.test import TestCase
from django.urls import reverse
from rest_framework import status
from rest_framework.test import APIClient
from core.models import (Reference, Species, Linelist,
                         SpeciesMetadata, MetaReference)
import tempfile
import os
import json
from rdkit import Chem
from rdkit.Chem import Descriptors
import selfies as sf


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


def create_species(**params):
    """Helper function to create a species."""
    defaults = {
        'name': ['common_name', 'some_name'],
        'iupac_name': 'Test IUPAC Name',
        'name_formula': 'Test Name Formula',
        'name_html': 'Test Name HTML',
        'molecular_mass': Descriptors.ExactMolWt(Chem.MolFromSmiles('C')),
        'smiles': 'C',
        'standard_inchi': 'Test InChI',
        'standard_inchi_key': 'Test InChI Key',
        'selfies': sf.encoder('C'),
        'mol_obj': 'C',
        'notes': 'Test Species',
    }
    defaults.update(params)
    defaults.update(molecular_mass=Descriptors.ExactMolWt(
        Chem.MolFromSmiles(defaults['smiles'])),
        selfies=sf.encoder(defaults['smiles']),
        mol_obj=defaults['smiles'])

    return Species.objects.create(**defaults)


def create_linelist(**params):
    """Helper function to create a linelist."""
    defaults = {
        'linelist_name': 'Test Linelist',
    }
    defaults.update(params)

    return Linelist.objects.create(**defaults)


def create_meta(species_id, linelist_id, **params):
    """Helper function to create a species metadata."""
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


class PublicReferenceApiTests(TestCase):
    """Test the publicly available reference API."""

    def setUp(self):
        self.client = APIClient()

    def test_auth_required_for_post(self):
        """Test that authentication is required for post
        creating references."""
        url = reverse('data:reference-list')
        payload = {
            'doi': '10.1021/acs.jcim.0c00128',
            'ref_url': 'https://doi.org/10.1021/acs.jcim.0c00128',
            'bibtex': 'bibtex_file',
            'notes': 'Test reference'}

        response = self.client.post(url, payload)
        self.assertEqual(response.status_code, status.HTTP_401_UNAUTHORIZED)

    def test_auth_required_for_put(self):
        """Test that authentication is required for put
        updating references."""
        reference = create_reference()
        url = reverse('data:reference-detail', args=[reference.id])
        payload = {
            'doi': '10.1021/acs.jcim.0c00128',
            'ref_url': 'https://doi.org/10.1021/acs.jcim.0c00128',
            'bibtex': 'bibtex_file',
            'notes': 'Test reference put',
            '_change_reason': 'Test change reason'}

        response = self.client.put(url, payload)
        self.assertEqual(response.status_code, status.HTTP_401_UNAUTHORIZED)

    def test_auth_required_for_patch(self):
        """Test that authentication is required for patch
        updating references."""
        reference = create_reference()
        url = reverse('data:reference-detail', args=[reference.id])
        payload = {'doi': 'new doi', '_change_reason': 'Test change reason'}

        response = self.client.patch(url, payload)
        self.assertEqual(response.status_code, status.HTTP_401_UNAUTHORIZED)

    def test_auth_required_for_delete(self):
        """Test that authentication is required for deleting references."""
        reference = create_reference()
        url = reverse('data:reference-detail',
                      args=[reference.id]) + \
            '?delete_reason=Test delete reason'

        response = self.client.delete(url)
        self.assertEqual(response.status_code, status.HTTP_401_UNAUTHORIZED)

    def test_get_reference_list(self):
        """Test getting a list of references."""
        create_reference(ref_url='url1')
        create_reference(ref_url='url2')
        url = reverse('data:reference-list')

        res = self.client.get(url)
        self.assertEqual(res.status_code, status.HTTP_200_OK)

    def test_get_reference_detail(self):
        """Test getting a reference detail."""
        reference = create_reference()
        url = reverse('data:reference-detail', args=[reference.id])

        res = self.client.get(url)
        self.assertEqual(res.status_code, status.HTTP_200_OK)


class PrivateReferenceApiTests(TestCase):
    """Test the private reference API."""

    def setUp(self):
        self.client = APIClient()
        self.user = get_user_model().objects.create_user(
            email='example@example.com',
            password='testpass',
            name='Test User',
            organization='Test Organization')
        self.client.force_authenticate(self.user)

    def test_create_reference(self):
        """Test creating a reference."""
        url = reverse('data:reference-list')
        with tempfile.NamedTemporaryFile(suffix='.bib') as bib_file:
            bib_file.write(b'@article{test, title={Test}}')
            bib_file.seek(0)
            payload = {
                'doi': 'test doi',
                'ref_url': 'test url',
                'bibtex': bib_file,
                'notes': 'Test reference'}
            res = self.client.post(url, payload, format='multipart')
        self.assertEqual(res.status_code, status.HTTP_201_CREATED)
        reference = Reference.objects.get(id=res.data['id'])
        for k, v in payload.items():
            if k != 'bibtex':
                self.assertEqual(v, getattr(reference, k))
        self.assertIn('bibtex', res.data)
        self.assertTrue(os.path.exists(reference.bibtex.path))
        history_exists = Reference.history.filter(
            id=reference.id).exists()
        self.assertTrue(history_exists)
        reference.bibtex.delete()

    def test_create_reference_with_duplicate_url_fails(self):
        """Test creating a reference with a duplicate ref_url fails."""
        create_reference(ref_url='test url')
        url = reverse('data:reference-list')
        with tempfile.NamedTemporaryFile(suffix='.bib') as bib_file:
            bib_file.write(b'@article{test, title={Test}}')
            bib_file.seek(0)
            payload = {
                'doi': 'test doi',
                'ref_url': 'test url',
                'bibtex': bib_file,
                'notes': 'Test reference'}
            res = self.client.post(url, payload, format='multipart')
        self.assertEqual(res.status_code, status.HTTP_400_BAD_REQUEST)

    def test_create_reference_with_invalid_bibtex(self):
        """Test creating a reference with an invalid bibtex file."""
        url = reverse('data:reference-list')
        payload = {
            'doi': 'test doi',
            'ref_url': 'test url',
            'bibtex': 'invalid bibtex file',
            'notes': 'Test reference'}
        res = self.client.post(url, payload, format='multipart')
        self.assertEqual(res.status_code, status.HTTP_400_BAD_REQUEST)

    def test_partial_update(self):
        """Test updating a reference with patch."""
        reference = create_reference(ref_url='url')
        url = reverse('data:reference-detail', args=[reference.id])
        payload = {'doi': 'new doi', '_change_reason': 'Test change reason'}
        res = self.client.patch(url, payload)
        self.assertEqual(res.status_code, status.HTTP_200_OK)
        reference.refresh_from_db()
        for k, v in payload.items():
            if k != '_change_reason':
                self.assertEqual(v, getattr(reference, k))
        self.assertEqual(Reference.history.filter(id=reference.id).first(
        ).history_change_reason, payload['_change_reason'])
        self.assertEqual(Reference.history.filter(id=reference.id).first(
        ).history_user_id, self.user.id)
        history_count = Reference.history.filter(id=reference.id).count()
        self.assertEqual(history_count, 2)

    def test_full_update(self):
        """Test updating a reference with put."""
        reference = create_reference(ref_url='Original reference')
        url = reverse('data:reference-detail', args=[reference.id])
        with tempfile.NamedTemporaryFile(suffix='.bib') as bib_file:
            bib_file.write(b'@article{test, title={Test}}')
            bib_file.seek(0)
            payload = {'doi': 'Updated doi',
                       'ref_url': 'Updated reference',
                       'bibtex': bib_file,
                       'notes': 'Updated notes',
                       '_change_reason': 'Test change reason'}
            res = self.client.put(url, payload, format='multipart')
        self.assertEqual(res.status_code, status.HTTP_200_OK)
        reference.refresh_from_db()
        for k, v in payload.items():
            if k not in ['_change_reason', 'bibtex']:
                self.assertEqual(v, getattr(reference, k))
        self.assertTrue(os.path.exists(reference.bibtex.path))
        self.assertEqual(Reference.history.filter(id=reference.id).first(
        ).history_change_reason, payload['_change_reason'])
        self.assertEqual(Reference.history.filter(id=reference.id).first(
        ).history_user_id, self.user.id)
        history_count = Reference.history.filter(id=reference.id).count()
        self.assertEqual(history_count, 2)
        reference.bibtex.delete()

    def test_full_update_ref_without_reason_fails(self):
        """Test updating a reference with put without change reason fails."""
        reference = create_reference(ref_url='Original reference')
        url = reverse('data:reference-detail', args=[reference.id])
        with tempfile.NamedTemporaryFile(suffix='.bib') as bib_file:
            bib_file.write(b'@article{test, title={Test}}')
            bib_file.seek(0)
            payload = {'doi': 'Updated doi',
                       'ref_url': 'Updated reference',
                       'bibtex': bib_file,
                       'notes': 'Updated notes'}
            res = self.client.put(url, payload, format='multipart')
        self.assertEqual(res.status_code, status.HTTP_400_BAD_REQUEST)

    def test_delete_reference(self):
        """Test deleting a reference."""
        reference = create_reference()
        url = reverse('data:reference-detail', args=[reference.id]) + \
            '?delete_reason=Test delete reason'
        res = self.client.delete(url)
        self.assertEqual(res.status_code, status.HTTP_204_NO_CONTENT)
        self.assertFalse(Reference.objects.filter(id=reference.id).exists())
        self.assertEqual(Reference.history.filter(id=reference.id).first(
        ).history_change_reason, 'Test delete reason')
        self.assertEqual(Reference.history.filter(id=reference.id).first(
        ).history_user_id, self.user.id)
        history_count = Reference.history.filter(id=reference.id).count()
        self.assertEqual(history_count, 2)

    def test_delete_reference_without_delete_reason_fails(self):
        """Test deleting a reference without delete reason fails."""
        reference = create_reference()
        url = reverse('data:reference-detail', args=[reference.id])
        res = self.client.delete(url)
        self.assertEqual(res.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertTrue(Reference.objects.filter(id=reference.id).exists())

    def test_delete_referenced_reference_fails(self):
        """Test deleting a referenced reference fails."""
        linelist = create_linelist()
        species = create_species()
        meta = create_meta(species.id, linelist.id)
        reference = create_reference()
        create_metaref(meta.id, reference.id)
        url = reverse('data:reference-detail', args=[reference.id]) + \
            '?delete_reason=Test delete reason'
        res = self.client.delete(url)
        self.assertEqual(res.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertTrue(Reference.objects.filter(id=reference.id).exists())

    def test_bibtex_without_id_fails(self):
        """Test getting merged bibtex files without id fails."""
        url = reverse('data:reference-get-bibtex')
        res = self.client.get(url)
        self.assertEqual(res.status_code, status.HTTP_400_BAD_REQUEST)
