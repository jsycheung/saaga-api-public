"""
Test species APIs.
"""
from django.contrib.auth import get_user_model
from django.test import TestCase
from django.urls import reverse
from rest_framework import status
from rest_framework.test import APIClient
from core.models import Species, SpeciesMetadata, Linelist
from rdkit import Chem
from rdkit.Chem import Descriptors
import selfies as sf
import json


def create_linelist(linelist_name='Test Linelist'):
    """Helper function to create a linelist."""
    return Linelist.objects.create(linelist_name=linelist_name)


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


class PublicSpeciesApiTests(TestCase):
    """Test the publicly available species API."""

    def setUp(self):
        self.client = APIClient()

    def test_auth_required_for_post(self):
        """Test that authentication is required for post creating species."""
        url = reverse('data:species-list')
        payload = {
            'name': ['common_name', 'test_name'],
            'iupac_name': 'Test IUPAC Name',
            'name_formula': 'Test Name Formula',
            'name_html': 'Test Name HTML',
            'smiles': 'C',
            'standard_inchi': 'Test InChI',
            'standard_inchi_key': 'Test InChI Key',
            'notes': 'Test Species',
        }

        response = self.client.post(url, payload)
        self.assertEqual(response.status_code, status.HTTP_401_UNAUTHORIZED)

    def test_auth_required_for_put(self):
        """Test that authentication is required for put updating species."""
        species = create_species()
        url = reverse('data:species-detail', args=[species.id])
        payload = {
            'name': ['common_name', 'test_name'],
            'iupac_name': 'Test IUPAC Name',
            'name_formula': 'Test Name Formula',
            'name_html': 'Test Name HTML',
            'smiles': 'C',
            'standard_inchi': 'Test InChI',
            'standard_inchi_key': 'Test InChI Key',
            'notes': 'Test Species',
            '_change_reason': 'Test change reason',
        }

        response = self.client.put(url, payload)
        self.assertEqual(response.status_code, status.HTTP_401_UNAUTHORIZED)

    def test_auth_required_for_patch(self):
        """Test that authentication is required for patch updating species."""
        species = create_species()
        url = reverse('data:species-detail', args=[species.id])
        payload = {'name': ['common_name1', 'test_name1'],
                   '_change_reason': 'Test change reason'}

        response = self.client.patch(url, payload)
        self.assertEqual(response.status_code, status.HTTP_401_UNAUTHORIZED)

    def test_auth_required_for_delete(self):
        """Test that authentication is required for deleting species."""
        species = create_species()
        url = reverse('data:species-detail',
                      args=[species.id]) + '?delete_reason=Test delete reason'

        response = self.client.delete(url)
        self.assertEqual(response.status_code, status.HTTP_401_UNAUTHORIZED)

    def test_get_species_list(self):
        """Test getting a list of species."""
        create_species(smiles='C', iupac_name='test iupac1')
        create_species(smiles='CC', iupac_name='test iupac2')
        url = reverse('data:species-list')

        res = self.client.get(url)
        self.assertEqual(res.status_code, status.HTTP_200_OK)
        self.assertEqual(len(res.data), 2)

    def test_get_species_detail(self):
        """Test getting a species detail."""
        species = create_species(smiles='CCC', iupac_name='test iupac3')
        url = reverse('data:species-detail', args=[species.id])

        res = self.client.get(url)
        self.assertEqual(res.status_code, status.HTTP_200_OK)


class PrivateSpeciesApiTests(TestCase):
    """Test the private species API."""

    def setUp(self):
        self.client = APIClient()
        self.user = get_user_model().objects.create_user(
            email='example@example.com',
            password='testpass',
            name='Test User',
            organization='Test Organization')
        self.client.force_authenticate(self.user)

    def test_create_species(self):
        """Test creating a species."""
        url = reverse('data:species-list')
        payload = {
            'name': json.dumps(["common_name", "Test Species"]),
            'iupac_name': 'Test IUPAC Name',
            'name_formula': 'Test Name Formula',
            'name_html': 'Test Name HTML',
            'smiles': 'CCCC',
            'standard_inchi': 'test inchi4',
            'standard_inchi_key': 'test inchi4',
            'notes': 'Test Species',
        }

        res = self.client.post(url, payload)
        self.assertEqual(res.status_code, status.HTTP_201_CREATED)
        species = Species.objects.get(id=res.data['id'])
        history_count = Species.history.filter(
            id=species.id).count()
        self.assertEqual(history_count, 1)

    def test_create_species_with_existing_iupac(self):
        """Test creating a species with an existing iupac name fails."""
        create_species(iupac_name='Test IUPAC Name existing')
        url = reverse('data:species-list')
        payload = {
            'name': json.dumps(['common_name', 'Test Species']),
            'iupac_name': 'Test IUPAC Name existing',
            'name_formula': 'Test Name Formula',
            'name_html': 'Test Name HTML',
            'smiles': 'CCCCC',
            'standard_inchi': 'test inchi5',
            'standard_inchi_key': 'CCCCC',
            'notes': 'Test Species',
        }

        res = self.client.post(url, payload)
        self.assertEqual(res.status_code, status.HTTP_400_BAD_REQUEST)

    def test_partial_update_species(self):
        """Test updating a species with patch."""
        species = create_species(iupac_name='test iupac partial update')
        url = reverse('data:species-detail', args=[species.id])
        payload = {'iupac_name': 'test name partial update',
                   '_change_reason': 'Test change reason'}

        res = self.client.patch(url, payload)
        self.assertEqual(res.status_code, status.HTTP_200_OK)
        species.refresh_from_db()
        self.assertEqual(species.iupac_name, payload['iupac_name'])
        self.assertEqual(Species.history.filter(id=species.id).first(
        ).history_change_reason, payload['_change_reason'])
        self.assertEqual(Species.history.filter(id=species.id).first(
        ).history_user_id, self.user.id)
        history_count = Species.history.filter(
            id=species.id).count()
        self.assertEqual(history_count, 2)

    def test_full_update_species(self):
        """Test updating a species with put."""
        species = create_species(iupac_name='test iupac full update')
        url = reverse('data:species-detail', args=[species.id])
        payload = {
            'name': json.dumps(['common_name', 'Test Species']),
            'iupac_name': 'Test IUPAC Name full update',
            'name_formula': 'Test Name Formula',
            'name_html': 'Test Name HTML',
            'smiles': 'CCCCC',
            'standard_inchi': 'test inchi full update payload',
            'standard_inchi_key': 'CCCCC',
            'molecular_mass': Descriptors.ExactMolWt(
                Chem.MolFromSmiles('CCCCC')),
            'selfies': sf.encoder('CCCCC'),
            'mol_obj': 'CCCCC',
            'notes': 'Test Species',
            '_change_reason': 'Test change reason',
        }

        res = self.client.put(url, payload)
        self.assertEqual(res.status_code, status.HTTP_200_OK)
        species.refresh_from_db()
        self.assertEqual(species.iupac_name, payload['iupac_name'])
        self.assertEqual(Species.history.filter(id=species.id).first(
        ).history_change_reason, payload['_change_reason'])
        self.assertEqual(Species.history.filter(id=species.id).first(
        ).history_user_id, self.user.id)
        history_count = Species.history.filter(
            id=species.id).count()
        self.assertEqual(history_count, 2)

    def test_full_update_species_without_reason_fails(self):
        """Test updating a species with put without change reason fails."""
        species = create_species(iupac_name='test iupac full update')
        url = reverse('data:species-detail', args=[species.id])
        payload = {
            'name': json.dumps(['common_name', 'Test Species']),
            'iupac_name': 'Test IUPAC Name full update',
            'name_formula': 'Test Name Formula',
            'name_html': 'Test Name HTML',
            'smiles': 'CCCCC',
            'standard_inchi': 'test inchi full update payload',
            'standard_inchi_key': 'CCCCC',
            'molecular_mass': Descriptors.ExactMolWt(
                Chem.MolFromSmiles('CCCCC')),
            'selfies': sf.encoder('CCCCC'),
            'mol_obj': 'CCCCC',
            'notes': 'Test Species',
        }

        res = self.client.put(url, payload)
        self.assertEqual(res.status_code, status.HTTP_400_BAD_REQUEST)

    def test_delete_species(self):
        """Test deleting a species."""
        species = create_species(iupac_name='test iupac delete')
        url = reverse('data:species-detail',
                      args=[species.id]) + \
            '?delete_reason=Test delete reason'

        res = self.client.delete(url)
        self.assertEqual(res.status_code, status.HTTP_204_NO_CONTENT)
        self.assertFalse(Species.objects.filter(id=species.id).exists())
        self.assertEqual(Species.history.filter(id=species.id).first(
        ).history_change_reason, 'Test delete reason')
        self.assertEqual(Species.history.filter(id=species.id).first(
        ).history_user_id, self.user.id)
        history_count = Species.history.filter(
            id=species.id).count()
        self.assertEqual(history_count, 2)

    def test_delete_species_without_delete_reason_fails(self):
        """Test deleting a species without delete reason fails."""
        species = create_species()
        url = reverse('data:species-detail', args=[species.id])
        res = self.client.delete(url)
        self.assertEqual(res.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertTrue(Species.objects.filter(id=species.id).exists())

    def test_delete_referenced_species_fails(self):
        """Test deleting a referenced species fails."""
        species = create_species()
        linelist = create_linelist()
        create_meta(species.id, linelist.id)
        url = reverse('data:species-detail',
                      args=[species.id]) + \
            '?delete_reason=Test delete reason'

        res = self.client.delete(url)
        self.assertEqual(res.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertTrue(Species.objects.filter(id=species.id).exists())
