"""
Test for Line APIs.
"""
from django.contrib.auth import get_user_model
from django.test import TestCase
from django.urls import reverse
from rest_framework import status
from rest_framework.test import APIClient
from core.models import Linelist, Species, SpeciesMetadata, Line
from data.serializers import LineSerializer
import json
from rdkit import Chem
from rdkit.Chem import Descriptors
import selfies as sf
import tempfile
from django.core.files import File as DjangoFile


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


def create_line(meta_id, **params):
    """Helper function to create a line."""
    defaults = {
        'meta_id': meta_id,
        'measured': False,
        'frequency': 100.000,
        'uncertainty': 0.001,
        'intensity': 0.001,
        's_ij_mu2': 1.0,
        'a_ij': 0.001,
        'lower_state_energy': 0.001,
        'upper_state_energy': 0.001,
        'lower_state_degeneracy': 1,
        'upper_state_degeneracy': 1,
        'lower_state_qn': json.dumps({'J': 1, 'Ka': 0, 'Kc': 0}),
        'upper_state_qn': json.dumps({'J': 1, 'Ka': 0, 'Kc': 1}),
        'rovibrational': False,
        'vib_qn': '',
        'pickett_qn_code': 303,
        'pickett_lower_state_qn': '010000',
        'pickett_upper_state_qn': '010001',
        'notes': 'test create line'
    }
    defaults.update(params)

    return Line.objects.create(**defaults)


class PublicLineApiTests(TestCase):
    """Test the publicly available line API."""

    def setUp(self):
        self.client = APIClient()

    def test_auth_required_for_post(self):
        """Test that authentication is required for post creating line."""
        url = reverse('data:line-list')
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        payload = {
            'meta': meta.id,
            'cat_file': 'test cat file',
            'qn_label_str': 'J,Ka,Kc',
            'contains_rovibrational': False,
            'vib_qn': '',
            'notes': 'Test post'}

        res = self.client.post(url, payload)
        self.assertEqual(res.status_code, status.HTTP_401_UNAUTHORIZED)

    def test_auth_required_for_put(self):
        """Test that authentication is required for put
        updating metareference."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        line = create_line(meta.id)
        url = reverse('data:line-detail', args=[line.id])
        payload = {
            'meta': meta.id,
            'cat_file': 'cat_file',
            'qn_label_str': 'J,Ka,Kc',
            'contains_rovibrational': False,
            'vib_qn': '',
            'notes': 'Test post',
            '_change_reason': 'Test change reason'
        }
        res = self.client.put(url, payload)
        self.assertEqual(res.status_code, status.HTTP_401_UNAUTHORIZED)

    def test_auth_required_for_patch(self):
        """Test that authentication is required for patch updating line."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        line = create_line(meta.id)
        url = reverse('data:line-detail', args=[line.id])
        payload = {'notes': 'test patch line',
                   '_change_reason': 'Test change reason'}
        res = self.client.patch(url, payload)
        self.assertEqual(res.status_code, status.HTTP_401_UNAUTHORIZED)

    def test_auth_required_for_delete(self):
        """Test that authentication is required for deleting line."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        line = create_line(meta.id)
        url = reverse('data:line-detail',
                      args=[line.id]) + '?delete_reason=Test delete reason'
        res = self.client.delete(url)
        self.assertEqual(res.status_code, status.HTTP_401_UNAUTHORIZED)

    def test_get_line_list(self):
        """Test getting a list of lines."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        create_line(meta.id, frequency=100.000)
        create_line(meta.id, frequency=200.000)
        url = reverse('data:line-list')

        res = self.client.get(url)
        lines = Line.objects.all().order_by('-id')
        serializer = LineSerializer(lines, many=True)
        self.assertEqual(res.status_code, status.HTTP_200_OK)
        self.assertEqual(res.data, serializer.data)

    def test_get_line_detail(self):
        """Test getting a line detail."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        create_line(meta.id)
        line = Line.objects.all().first()
        url = reverse('data:line-detail', args=[line.id])

        res = self.client.get(url)
        serializer = LineSerializer(line)
        self.assertEqual(res.status_code, status.HTTP_200_OK)
        self.assertEqual(res.data, serializer.data)


class PrivateLineApiTests(TestCase):
    """Test the private line API."""

    def setUp(self):
        self.client = APIClient()
        self.user = get_user_model().objects.create_user(
            email='example@example.com',
            password='testpass',
            name='Test User',
            organization='Test Organization')
        self.client.force_authenticate(self.user)

    def test_create_line_with_wrong_cat_file_fails(self):
        """Test creating a line with wrong cat file fails."""
        url = reverse('data:line-list')
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        not_cat_file = tempfile.NamedTemporaryFile(suffix='.notcat')
        payload = {
            'meta': meta.id,
            'cat_file': DjangoFile(not_cat_file),
            'qn_label_str': 'J,Ka,Kc',
            'contains_rovibrational': False,
            'vib_qn': '',
            'notes': 'Test post'}

        res = self.client.post(url, payload)
        self.assertEqual(res.status_code, status.HTTP_400_BAD_REQUEST)

    def test_create_line_with_rovibrational_no_vib_qn_fails(self):
        """Test creating a line with rovibrational true but no vib qn fails"""
        url = reverse('data:line-list')
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        cat_file = tempfile.NamedTemporaryFile(suffix='.cat')
        payload = {
            'meta': meta.id,
            'cat_file': DjangoFile(cat_file),
            'qn_label_str': 'J,Ka,Kc',
            'contains_rovibrational': True,
            'vib_qn': '',
            'notes': 'Test post'}

        res = self.client.post(url, payload)
        self.assertEqual(res.status_code, status.HTTP_400_BAD_REQUEST)

    def test_create_line_with_rovibrational_wrong_vib_qn_fails(self):
        """Test creating a line with rovibrational wrong but vib qn fails."""
        url = reverse('data:line-list')
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        cat_file = tempfile.NamedTemporaryFile(suffix='.cat')
        payload = {
            'meta': meta.id,
            'cat_file': DjangoFile(cat_file),
            'qn_label_str': 'J,Ka,Kc',
            'contains_rovibrational': True,
            'vib_qn': 'N',
            'notes': 'Test post'}

        res = self.client.post(url, payload)
        self.assertEqual(res.status_code, status.HTTP_400_BAD_REQUEST)

    def test_create_line_not_rovibrational_with_vib_qn_fails(self):
        """Test creating line that is not rovibrational but with
        vib qn provided fails."""
        url = reverse('data:line-list')
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        cat_file = tempfile.NamedTemporaryFile(suffix='.cat')
        payload = {
            'meta': meta.id,
            'cat_file': DjangoFile(cat_file),
            'qn_label_str': 'J,Ka,Kc',
            'contains_rovibrational': False,
            'vib_qn': 'Kc',
            'notes': 'Test post'}

        res = self.client.post(url, payload)
        self.assertEqual(res.status_code, status.HTTP_400_BAD_REQUEST)

    def test_partial_update(self):
        """Test updating a line with patch."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        line = create_line(meta.id)
        url = reverse('data:line-detail', args=[line.id])
        payload = {'notes': 'Updated line patch',
                   '_change_reason': 'Test patch'}
        res = self.client.patch(url, payload)
        self.assertEqual(res.status_code, status.HTTP_200_OK)
        line.refresh_from_db()
        self.assertEqual(Line.history.filter(id=line.id).first(
        ).history_change_reason, payload['_change_reason'])
        self.assertEqual(Line.history.filter(id=line.id).first(
        ).history_user_id, self.user.id)
        history_count = Line.history.filter(id=line.id).count()
        self.assertEqual(history_count, 2)

    def test_full_update(self):
        """Test updating a line with put."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        line = create_line(meta.id)
        url = reverse('data:line-detail', args=[line.id])
        payload = {
            'meta': meta.id,
            'measured': False,
            'frequency': 100.000,
            'uncertainty': 0.001,
            'intensity': 0.001,
            's_ij_mu2': 1.0,
            'a_ij': 0.001,
            'lower_state_energy': 0.001,
            'upper_state_energy': 0.001,
            'lower_state_degeneracy': 1,
            'upper_state_degeneracy': 1,
            'lower_state_qn': json.dumps({'J': 1, 'Ka': 0, 'Kc': 0}),
            'upper_state_qn': json.dumps({'J': 1, 'Ka': 0, 'Kc': 1}),
            'rovibrational': False,
            'vib_qn': '',
            'pickett_qn_code': 303,
            'pickett_lower_state_qn': '010000',
            'pickett_upper_state_qn': '010001',
            'notes': 'test create line',
            '_change_reason': 'Test put full update'
        }
        res = self.client.put(url, payload)
        self.assertEqual(res.status_code, status.HTTP_200_OK)
        line.refresh_from_db()
        self.assertEqual(Line.history.filter(id=line.id).first(
        ).history_change_reason, payload['_change_reason'])
        self.assertEqual(Line.history.filter(id=line.id).first(
        ).history_user_id, self.user.id)
        history_count = Line.history.filter(id=line.id).count()
        self.assertEqual(history_count, 2)

    def test_full_update_line_without_reason_fails(self):
        """Test updating a line with put without change reason fails."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        line = create_line(meta.id)
        url = reverse('data:line-detail', args=[line.id])
        payload = {
            'meta': meta.id,
            'measured': False,
            'frequency': 100.000,
            'uncertainty': 0.001,
            'intensity': 0.001,
            's_ij_mu2': 1.0,
            'a_ij': 0.001,
            'lower_state_energy': 0.001,
            'upper_state_energy': 0.001,
            'lower_state_degeneracy': 1,
            'upper_state_degeneracy': 1,
            'lower_state_qn': json.dumps({'J': 1, 'Ka': 0, 'Kc': 0}),
            'upper_state_qn': json.dumps({'J': 1, 'Ka': 0, 'Kc': 1}),
            'rovibrational': False,
            'vib_qn': '',
            'pickett_qn_code': 303,
            'pickett_lower_state_qn': '010000',
            'pickett_upper_state_qn': '010001',
            'notes': 'test create line',
        }
        res = self.client.put(url, payload)
        self.assertEqual(res.status_code, status.HTTP_400_BAD_REQUEST)

    def test_delete_line(self):
        """Test deleting a line."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        line = create_line(meta.id)
        url = reverse('data:line-detail',
                      args=[line.id]) + \
            '?delete_reason=Test delete reason line'
        res = self.client.delete(url)
        self.assertEqual(res.status_code, status.HTTP_204_NO_CONTENT)
        self.assertFalse(Line.objects.filter(id=line.id).exists())
        self.assertEqual(Line.history.filter(id=line.id).first(
        ).history_change_reason, 'Test delete reason line')
        self.assertEqual(Line.history.filter(id=line.id).first(
        ).history_user_id, self.user.id)
        history_count = Line.history.filter(id=line.id).count()
        self.assertEqual(history_count, 2)

    def delete_line_without_delete_reason_fails(self):
        """Test deleting a line without delete reason fails."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        line = create_line(meta.id)
        url = reverse('data:line-detail', args=[line.id])
        res = self.client.delete(url)
        self.assertEqual(res.status_code, status.HTTP_400_BAD_REQUEST)

    def test_query_lines_min_max(self):
        """Test querying lines with min and max frequencies both specified."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        create_line(meta.id, frequency=100.000)
        create_line(meta.id, frequency=200.000)
        url = reverse('data:line-query') + '?min_freq=99.000&max_freq=201.000'
        res = self.client.get(url)
        self.assertEqual(res.status_code, status.HTTP_200_OK)
        self.assertEqual(len(res.data), 2)

    def test_query_lines_min(self):
        """Test querying lines with min frequency specified."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        create_line(meta.id, frequency=100.000)
        create_line(meta.id, frequency=200.000)
        url = reverse('data:line-query') + '?min_freq=99.000'
        res = self.client.get(url)
        self.assertEqual(res.status_code, status.HTTP_200_OK)
        self.assertEqual(len(res.data), 2)

    def test_query_lines_max(self):
        """Test querying lines with max frequency specified."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        create_line(meta.id, frequency=100.000)
        create_line(meta.id, frequency=200.000)
        url = reverse('data:line-query') + '?max_freq=101.000'
        res = self.client.get(url)
        self.assertEqual(res.status_code, status.HTTP_200_OK)
        self.assertEqual(len(res.data), 1)

    def test_query_lines_without_freq_fails(self):
        """Test querying lines without frequency specified fails."""
        url = reverse('data:line-query')
        res = self.client.get(url)
        self.assertEqual(res.status_code, status.HTTP_400_BAD_REQUEST)
