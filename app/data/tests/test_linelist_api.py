from django.contrib.auth import get_user_model
from django.test import TestCase
from django.urls import reverse
from rest_framework import status
from rest_framework.test import APIClient
from core.models import Linelist
from data.serializers import LinelistSerializer, Species, SpeciesMetadata
from rdkit import Chem
from rdkit.Chem import Descriptors
import json
import selfies as sf
from django.db.utils import IntegrityError


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


class PublicLinelistApiTests(TestCase):
    """Test the publicly available linelist API."""

    def setUp(self):
        self.client = APIClient()

    def test_auth_required_for_post(self):
        """Test that authentication is required for post
        creating linelists."""
        url = reverse('data:linelist-list')
        payload = {'linelist_name': 'Test Linelist'}

        response = self.client.post(url, payload)
        self.assertEqual(response.status_code, status.HTTP_401_UNAUTHORIZED)

    def test_auth_required_for_put(self):
        """Test that authentication is required for put updating linelists."""
        linelist = create_linelist()
        url = reverse('data:linelist-detail', args=[linelist.id])
        payload = {'linelist_name': 'Updated Linelist'}

        response = self.client.put(url, payload)
        self.assertEqual(response.status_code, status.HTTP_401_UNAUTHORIZED)

    def test_auth_required_for_patch(self):
        """Test that authentication is required for
        patch updating linelists."""
        linelist = create_linelist()
        url = reverse('data:linelist-detail', args=[linelist.id])
        payload = {'linelist_name': 'Updated Linelist'}

        response = self.client.patch(url, payload)
        self.assertEqual(response.status_code, status.HTTP_401_UNAUTHORIZED)

    def test_auth_required_for_delete(self):
        """Test that authentication is required for deleting linelists."""
        linelist = create_linelist()
        url = reverse('data:linelist-detail', args=[linelist.id])

        response = self.client.delete(url)
        self.assertEqual(response.status_code, status.HTTP_401_UNAUTHORIZED)

    def test_get_linelist_list(self):
        """Test getting a list of linelists."""
        create_linelist(linelist_name='Test Linelist 1')
        create_linelist(linelist_name='Test Linelist 2')
        url = reverse('data:linelist-list')

        res = self.client.get(url)
        linelists = Linelist.objects.all().order_by('-id')
        serializer = LinelistSerializer(linelists, many=True)
        self.assertEqual(res.status_code, status.HTTP_200_OK)
        self.assertEqual(res.data, serializer.data)

    def test_get_linelist_detail(self):
        """Test getting a linelist detail."""
        linelist = create_linelist()
        url = reverse('data:linelist-detail', args=[linelist.id])

        res = self.client.get(url)
        serializer = LinelistSerializer(linelist)
        self.assertEqual(res.status_code, status.HTTP_200_OK)
        self.assertEqual(res.data, serializer.data)


class PrivateLinelistApiTests(TestCase):
    """Test the private linelist API."""

    def setUp(self):
        self.client = APIClient()
        self.user = get_user_model().objects.create_user(
            email='example@example.com',
            password='testpass',
            name='Test User',
            organization='Test Organization')
        self.client.force_authenticate(self.user)

    def test_create_linelist(self):
        """Test creating a linelist."""
        payload = {'linelist_name': 'Test Linelist'}
        url = reverse('data:linelist-list')
        res = self.client.post(url, payload)
        self.assertEqual(res.status_code, status.HTTP_201_CREATED)
        linelist = Linelist.objects.get(id=res.data['id'])
        for k, v in payload.items():
            self.assertEqual(v.lower(), getattr(linelist, k))
        history_exists = Linelist.history.filter(
            linelist_name=payload['linelist_name'].lower()).exists()
        self.assertTrue(history_exists)

    def test_create_duplicate_linelist_fails(self):
        """Test creating a duplicate linelist fails."""
        create_linelist(linelist_name='Test Linelist duplicate')
        payload = {'linelist_name': 'Test Linelist duplicate'}
        url = reverse('data:linelist-list')
        with self.assertRaises(IntegrityError):
            self.client.post(url, payload)

    def test_partial_update(self):
        """Test updating a linelist with patch."""
        linelist = create_linelist(linelist_name='Original linelist')
        url = reverse('data:linelist-detail', args=[linelist.id])
        payload = {'linelist_name': 'Updated Linelist',
                   '_change_reason': 'Test patch'}
        res = self.client.patch(url, payload)
        self.assertEqual(res.status_code, status.HTTP_200_OK)
        linelist.refresh_from_db()
        for k, v in payload.items():
            if k != '_change_reason':
                self.assertEqual(v.lower(), getattr(linelist, k))
        self.assertEqual(Linelist.history.filter(id=linelist.id).first(
        ).history_change_reason, payload['_change_reason'])
        self.assertEqual(Linelist.history.filter(id=linelist.id).first(
        ).history_user_id, self.user.id)
        history_count = Linelist.history.filter(id=linelist.id).count()
        self.assertEqual(history_count, 2)

    def test_full_update(self):
        """Test updating a linelist with put."""
        linelist = create_linelist(
            linelist_name='Original linelist')
        url = reverse('data:linelist-detail', args=[linelist.id])
        payload = {'linelist_name': 'Updated Linelist',
                   '_change_reason': 'Test put'}
        res = self.client.put(url, payload)
        self.assertEqual(res.status_code, status.HTTP_200_OK)
        linelist.refresh_from_db()
        for k, v in payload.items():
            if k != '_change_reason':
                self.assertEqual(v.lower(), getattr(linelist, k))
        self.assertEqual(Linelist.history.filter(id=linelist.id).first(
        ).history_change_reason, payload['_change_reason'])
        self.assertEqual(Linelist.history.filter(id=linelist.id).first(
        ).history_user_id, self.user.id)
        history_count = Linelist.history.filter(id=linelist.id).count()
        self.assertEqual(history_count, 2)

    def test_full_update_linelist_without_reason_fails(self):
        """Test updating a linelist with put without a change reason fails."""
        linelist = create_linelist(
            linelist_name='Original linelist')
        url = reverse('data:linelist-detail', args=[linelist.id])
        payload = {'linelist_name': 'Updated Linelist'}
        res = self.client.put(url, payload)
        self.assertEqual(res.status_code, status.HTTP_400_BAD_REQUEST)

    def test_delete_linelist(self):
        """Test deleting a linelist."""
        linelist = create_linelist()
        url = reverse('data:linelist-detail',
                      args=[linelist.id]) + '?delete_reason=Test delete'
        res = self.client.delete(url)
        self.assertEqual(res.status_code, status.HTTP_204_NO_CONTENT)
        self.assertFalse(Linelist.objects.filter(id=linelist.id).exists())
        self.assertEqual(Linelist.history.filter(id=linelist.id).first(
        ).history_change_reason, 'Test delete')
        self.assertEqual(Linelist.history.filter(id=linelist.id).first(
        ).history_user_id, self.user.id)
        history_count = Linelist.history.filter(id=linelist.id).count()
        self.assertEqual(history_count, 2)

    def test_delete_linelist_without_reason_fails(self):
        """Test deleting a linelist without a delete reason fails."""
        linelist = create_linelist()
        url = reverse('data:linelist-detail', args=[linelist.id])
        res = self.client.delete(url)
        self.assertEqual(res.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertTrue(Linelist.objects.filter(id=linelist.id).exists())

    def test_delete_referenced_linelist_fails(self):
        """Test deleting a linelist that is referenced by species
        metadata fails."""
        linelist = create_linelist()
        species = create_species()
        create_meta(species.id, linelist.id)
        url = reverse('data:linelist-detail',
                      args=[linelist.id]) + \
            '?delete_reason=Test delete referenced'
        res = self.client.delete(url)
        self.assertEqual(res.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertTrue(Linelist.objects.filter(id=linelist.id).exists())
