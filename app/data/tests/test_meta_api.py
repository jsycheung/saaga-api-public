"""
Test species metadata APIs.
"""
from django.contrib.auth import get_user_model
from django.test import TestCase
from django.urls import reverse
from rest_framework import status
from rest_framework.test import APIClient
from core.models import SpeciesMetadata, Species, Linelist, Line
import json
from rdkit import Chem
from rdkit.Chem import Descriptors
import selfies as sf
import tempfile


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


class PublicMetaApiTests(TestCase):
    """Test the publicly available species metadata API."""

    def setUp(self):
        self.client = APIClient()

    def test_auth_required_for_post(self):
        """Test that authentication is required for creating
        species metadata."""
        url = reverse('data:speciesmetadata-list')
        species = create_species()
        linelist = create_linelist()
        payload = {
            'species_id': species.id,
            'molecule_tag': 1,
            'hyperfine': False,
            'degree_of_freedom': 3,
            'category': 'asymmetric top',
            'mu_a': 0.5,
            'mu_b': 0,
            'mu_c': 0,
            'a_const': 0.123,
            'b_const': 0,
            'c_const': 0,
            'linelist_id': linelist.id,
            'data_date': '2020-01-01',
            'data_contributor': 'Test Contributor',
            'int_file': 'test_int_file',
            'var_file': 'test_var_file',
            'fit_file': 'test_fit_file',
            'lin_file': 'test_lin_file',
            'qpart_file': 'test_qpart_file',
            'notes': 'Test Species Metadata',
        }

        response = self.client.post(url, payload)
        self.assertEqual(response.status_code, status.HTTP_401_UNAUTHORIZED)

    def test_auth_required_for_put(self):
        """Test that authentication is required for put updating
        species metadata."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        url = reverse('data:speciesmetadata-detail', args=[meta.id])
        payload = {
            'species_id': species.id,
            'molecule_tag': 1,
            'hyperfine': True,
            'degree_of_freedom': 3,
            'category': 'asymmetric top',
            'mu_a': 0.5,
            'mu_b': 0,
            'mu_c': 0,
            'a_const': 0.123,
            'b_const': 0,
            'c_const': 0,
            'linelist_id': linelist.id,
            'data_date': '2020-01-01',
            'data_contributor': 'Test Contributor',
            'int_file': 'test_int_file',
            'var_file': 'test_var_file',
            'fit_file': 'test_fit_file',
            'lin_file': 'test_lin_file',
            'qpart_file': 'test_qpart_file',
            'notes': 'Test Species Metadata',
            '_change_reason': 'Test change reason',
        }

        response = self.client.put(url, payload)
        self.assertEqual(response.status_code, status.HTTP_401_UNAUTHORIZED)

    def test_auth_required_for_patch(self):
        """Test that authentication is required for patch updating
        species metadata."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        url = reverse('data:speciesmetadata-detail', args=[meta.id])
        payload = {'degree_of_freedom': 2,
                   '_change_reason': 'Test change reason patch'}

        response = self.client.patch(url, payload)
        self.assertEqual(response.status_code, status.HTTP_401_UNAUTHORIZED)

    def test_auth_required_for_delete(self):
        """Test that authentication is required for deleting
        species metadata."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        url = reverse('data:speciesmetadata-detail',
                      args=[meta.id]) + '?delete_reason=Test delete reason'

        response = self.client.delete(url)
        self.assertEqual(response.status_code, status.HTTP_401_UNAUTHORIZED)

    def test_get_meta_list(self):
        """Test getting a list of species metadata."""
        species = create_species()
        linelist1 = create_linelist(linelist_name='Test Linelist 1')
        linelist2 = create_linelist(linelist_name='Test Linelist 2')
        create_meta(species.id, linelist1.id)
        create_meta(species.id, linelist2.id)
        url = reverse('data:speciesmetadata-list')

        res = self.client.get(url)
        self.assertEqual(res.status_code, status.HTTP_200_OK)
        self.assertEqual(len(res.data), 2)

    def test_get_meta_detail(self):
        """Test getting a species metadata detail."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        url = reverse('data:speciesmetadata-detail', args=[meta.id])

        res = self.client.get(url)
        self.assertEqual(res.status_code, status.HTTP_200_OK)


class PrivateSpeciesApiTests(TestCase):
    """Test the private species metadata API."""

    def setUp(self):
        self.client = APIClient()
        self.user = get_user_model().objects.create_user(
            email='example@example.com',
            password='testpass',
            name='Test User',
            organization='Test Organization')
        self.client.force_authenticate(self.user)

    def test_create_meta(self):
        """Test creating a species metadata."""
        species = create_species()
        linelist = create_linelist()
        url = reverse('data:speciesmetadata-list')
        with tempfile.NamedTemporaryFile(suffix='.int') as int_file, \
                tempfile.NamedTemporaryFile(suffix='.var') as var_file, \
                tempfile.NamedTemporaryFile(suffix='.fit') as fit_file, \
                tempfile.NamedTemporaryFile(suffix='.lin') as lin_file, \
                tempfile.NamedTemporaryFile(suffix='.qpart') as qpart_file:
            int_file.write(b"""FURAN
0  91 43282.2944  00  150  -7  -7 150.0  5
 001  0.661 /muA
 002  0.000 /muB
 003  0.000 /muC""")
            int_file.seek(0)
            var_file.write(b"FURAN                        "
                           b"                "
                           b"WedThu MaWeFri May 19 17:01:15 2023\n"
                           b"   3    9   40    0     0.0000E+00"
                           b"     5.0000E+20     "
                           b"1.0000E+00 1.0000000000\n"
                           b"S   1  1  0  50  0  1  1  1  1  -1   0\n"
                           b"        10000   9.447840918743423E+03 "
                           b"6.36184172E-02 /A\n"
                           b"        20000   9.246852060423453E+03 "
                           b"1.63997319E-02 /B\n"
                           b"        30000   4.671341872421471E+03 "
                           b"5.27198401E-02 /C")
            var_file.seek(0)
            fit_file.write(b'I am a fit file')
            fit_file.seek(0)
            lin_file.write(b'I am a lin file')
            lin_file.seek(0)
            qpart_file.write(b"""#form : interpolation
300.000   331777.6674""")
            qpart_file.seek(0)
            payload = {
                'species': species.id,
                'molecule_tag': 1,
                'hyperfine': True,
                'degree_of_freedom': 3,
                'category': 'asymmetric top',
                'mu_a': 0.5,
                'mu_b': 0,
                'mu_c': 0,
                'a_const': 0.123,
                'b_const': 0,
                'c_const': 0,
                'linelist': linelist.id,
                'data_date': '2020-01-01',
                'data_contributor': 'Test Contributor',
                'int_file': int_file,
                'var_file': var_file,
                'fit_file': fit_file,
                'lin_file': lin_file,
                'qpart_file': qpart_file,
                'notes': 'Test Species Metadata',
            }
            res = self.client.post(url, payload, format='multipart')

        self.assertEqual(res.status_code, status.HTTP_201_CREATED)
        meta = SpeciesMetadata.objects.get(id=res.data['id'])
        history_count = SpeciesMetadata.history.filter(
            id=meta.id).count()
        self.assertEqual(history_count, 1)
        meta.int_file.delete()
        meta.var_file.delete()
        meta.fit_file.delete()
        meta.lin_file.delete()
        meta.qpart_file.delete()

    def test_create_meta_without_300K(self):
        """Test creating a species metadata with qpart
        that does not contain 300K fails."""
        species = create_species()
        linelist = create_linelist()
        url = reverse('data:speciesmetadata-list')
        with tempfile.NamedTemporaryFile(suffix='.int') as int_file, \
                tempfile.NamedTemporaryFile(suffix='.var') as var_file, \
                tempfile.NamedTemporaryFile(suffix='.fit') as fit_file, \
                tempfile.NamedTemporaryFile(suffix='.lin') as lin_file, \
                tempfile.NamedTemporaryFile(suffix='.qpart') as qpart_file:
            int_file.write(b"FURAN"
                           b"0  91 43282.2944  00  150  -7  -7 150.0  5\n"
                           b" 001  0.661 /muA\n"
                           b" 002  0.000 /muB\n"
                           b" 003  0.000 /muC")
            int_file.seek(0)
            var_file.write(b"FURAN                                        "
                           b"WedThu MaWeFri May 19 17:01:15 2023\n"
                           b"   3    9   40    0     0.0000E+00"
                           b"     5.0000E+20     "
                           b"1.0000E+00 1.0000000000\n"
                           b"S   1  1  0  50  0  1  1  1  1  -1   0\n"
                           b"        10000   9.447840918743423E+03 "
                           b"6.36184172E-02 /A\n"
                           b"        20000   9.246852060423453E+03 "
                           b"1.63997319E-02 /B\n"
                           b"        30000   4.671341872421471E+03 "
                           b"5.27198401E-02 /C")
            var_file.seek(0)
            fit_file.write(b'I am a fit file')
            fit_file.seek(0)
            lin_file.write(b'I am a lin file')
            lin_file.seek(0)
            qpart_file.write(b"#form : interpolation\n"
                             b"200.000   331777.6674")
            qpart_file.seek(0)
            payload = {
                'species': species.id,
                'molecule_tag': 1,
                'hyperfine': True,
                'degree_of_freedom': 3,
                'category': 'asymmetric top',
                'mu_a': 0.5,
                'mu_b': 0,
                'mu_c': 0,
                'a_const': 0.123,
                'b_const': 0,
                'c_const': 0,
                'linelist': linelist.id,
                'data_date': '2020-01-01',
                'data_contributor': 'Test Contributor',
                'int_file': int_file,
                'var_file': var_file,
                'fit_file': fit_file,
                'lin_file': lin_file,
                'qpart_file': qpart_file,
                'notes': 'Test Species Metadata',
            }
            res = self.client.post(url, payload, format='multipart')

        self.assertEqual(res.status_code, status.HTTP_400_BAD_REQUEST)

    def test_create_meta_with_wrong_file_fails(self):
        """Test creating a species metadata with wrong file
        extension fails."""
        species = create_species()
        linelist = create_linelist()
        url = reverse('data:speciesmetadata-list')
        with tempfile.NamedTemporaryFile(suffix='.notint') as not_int_file, \
                tempfile.NamedTemporaryFile(suffix='.var') as var_file, \
                tempfile.NamedTemporaryFile(suffix='.fit') as fit_file, \
                tempfile.NamedTemporaryFile(suffix='.lin') as lin_file, \
                tempfile.NamedTemporaryFile(suffix='.qpart') as qpart_file:
            not_int_file.write(b"FURAN\n"
                               b"0  91 43282.2944  00  150  -7  -7 150.0  5\n"
                               b" 001  0.661 /muA\n"
                               b" 002  0.000 /muB\n"
                               b" 003  0.000 /muC")
            not_int_file.seek(0)
            var_file.write(b"FURAN                                        "
                           b"WedThu MaWeFri May 19 17:01:15 2023\n"
                           b"   3    9   40    0     0.0000E+00     "
                           b"5.0000E+20     "
                           b"1.0000E+00 1.0000000000\n"
                           b"S   1  1  0  50  0  1  1  1  1  -1   0\n"
                           b"        10000   9.447840918743423E+03 "
                           b"6.36184172E-02 /A\n"
                           b"        20000   9.246852060423453E+03 "
                           b"1.63997319E-02 /B\n"
                           b"        30000   4.671341872421471E+03 "
                           b"5.27198401E-02 /C")
            var_file.seek(0)
            fit_file.write(b'I am a fit file')
            fit_file.seek(0)
            lin_file.write(b'I am a lin file')
            lin_file.seek(0)
            qpart_file.write(b"#form : interpolation\n"
                             b"200.000   331777.6674")
            qpart_file.seek(0)
            payload = {
                'species': species.id,
                'molecule_tag': 1,
                'hyperfine': True,
                'degree_of_freedom': 3,
                'category': 'asymmetric top',
                'mu_a': 0.5,
                'mu_b': 0,
                'mu_c': 0,
                'a_const': 0.123,
                'b_const': 0,
                'c_const': 0,
                'linelist': linelist.id,
                'data_date': '2020-01-01',
                'data_contributor': 'Test Contributor',
                'int_file': not_int_file,
                'var_file': var_file,
                'fit_file': fit_file,
                'lin_file': lin_file,
                'qpart_file': qpart_file,
                'notes': 'Test Species Metadata',
            }
            res = self.client.post(url, payload, format='multipart')

        self.assertEqual(res.status_code, status.HTTP_400_BAD_REQUEST)

    def test_partial_update_meta(self):
        """Test updating a species metadata with patch."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        url = reverse('data:speciesmetadata-detail', args=[meta.id])
        payload = {'degree_of_freedom': 2,
                   '_change_reason': 'Test change reason'}

        res = self.client.patch(url, payload)
        self.assertEqual(res.status_code, status.HTTP_200_OK)
        meta.refresh_from_db()
        self.assertEqual(meta.degree_of_freedom, payload['degree_of_freedom'])
        self.assertEqual(SpeciesMetadata.history.filter(id=meta.id).first(
        ).history_change_reason, payload['_change_reason'])
        self.assertEqual(SpeciesMetadata.history.filter(id=meta.id).first(
        ).history_user_id, self.user.id)
        history_count = SpeciesMetadata.history.filter(
            id=meta.id).count()
        self.assertEqual(history_count, 2)

    def test_full_update_meta(self):
        """Test updating a species metadata with put."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        url = reverse('data:speciesmetadata-detail', args=[meta.id])
        with tempfile.NamedTemporaryFile(suffix='.int') as int_file, \
                tempfile.NamedTemporaryFile(suffix='.var') as var_file, \
                tempfile.NamedTemporaryFile(suffix='.fit') as fit_file, \
                tempfile.NamedTemporaryFile(suffix='.lin') as lin_file, \
                tempfile.NamedTemporaryFile(suffix='.qpart') as qpart_file:
            int_file.write(b"FURAN\n"
                           b"0  91 43282.2944  00  150  -7  -7 150.0  5\n"
                           b" 001  0.661 /muA\n"
                           b" 002  0.000 /muB\n"
                           b" 003  0.000 /muC")
            int_file.seek(0)
            var_file.write(b"FURAN                                        "
                           b"WedThu MaWeFri May 19 17:01:15 2023\n"
                           b"   3    9   40    0     0.0000E+00     "
                           b"5.0000E+20     "
                           b"1.0000E+00 1.0000000000\n"
                           b"S   1  1  0  50  0  1  1  1  1  -1   0\n"
                           b"        10000   9.447840918743423E+03 "
                           b"6.36184172E-02 /A\n"
                           b"        20000   9.246852060423453E+03 "
                           b"1.63997319E-02 /B\n"
                           b"        30000   4.671341872421471E+03 "
                           b"5.27198401E-02 /C")
            var_file.seek(0)
            fit_file.write(b'I am a fit file')
            fit_file.seek(0)
            lin_file.write(b'I am a lin file')
            lin_file.seek(0)
            qpart_file.write(b"#form : interpolation\n"
                             b"200.000   331777.6674")
            qpart_file.seek(0)
            payload = {
                'species': species.id,
                'molecule_tag': 1,
                'hyperfine': True,
                'degree_of_freedom': 3,
                'category': 'asymmetric top',
                'partition_function': json.dumps({'300.000': '331777.6674'}),
                'mu_a': 0.5,
                'mu_b': 0,
                'mu_c': 0,
                'a_const': 0.123,
                'b_const': 0,
                'c_const': 0,
                'linelist': linelist.id,
                'data_date': '2020-01-01',
                'data_contributor': 'Test Contributor',
                'int_file': int_file,
                'var_file': var_file,
                'fit_file': fit_file,
                'lin_file': lin_file,
                'qpart_file': qpart_file,
                'notes': 'Test Species Metadata Put',
                '_change_reason': 'Test change reason',
            }
            res = self.client.put(url, payload, format='multipart')

        self.assertEqual(res.status_code, status.HTTP_200_OK)
        meta.refresh_from_db()
        self.assertEqual(meta.notes, payload['notes'])
        self.assertEqual(SpeciesMetadata.history.filter(id=meta.id).first(
        ).history_change_reason, payload['_change_reason'])
        self.assertEqual(SpeciesMetadata.history.filter(id=meta.id).first(
        ).history_user_id, self.user.id)
        history_count = SpeciesMetadata.history.filter(
            id=meta.id).count()
        self.assertEqual(history_count, 2)
        meta.int_file.delete()
        meta.var_file.delete()
        meta.fit_file.delete()
        meta.lin_file.delete()
        meta.qpart_file.delete()

    def test_full_update_meta_without_reason_fails(self):
        """Test updating a species metadata with put without
        change reason fails."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        url = reverse('data:speciesmetadata-detail', args=[meta.id])
        with tempfile.NamedTemporaryFile(suffix='.int') as int_file, \
                tempfile.NamedTemporaryFile(suffix='.var') as var_file, \
                tempfile.NamedTemporaryFile(suffix='.fit') as fit_file, \
                tempfile.NamedTemporaryFile(suffix='.lin') as lin_file, \
                tempfile.NamedTemporaryFile(suffix='.qpart') as qpart_file:
            int_file.write(b"FURAN\n"
                           b"0  91 43282.2944  00  150  -7  -7 150.0  5\n"
                           b" 001  0.661 /muA\n"
                           b" 002  0.000 /muB\n"
                           b" 003  0.000 /muC")
            int_file.seek(0)
            var_file.write(b"FURAN                                        "
                           b"WedThu MaWeFri May 19 17:01:15 2023\n"
                           b"   3    9   40    0     0.0000E+00     "
                           b"5.0000E+20     "
                           b"1.0000E+00 1.0000000000\n"
                           b"S   1  1  0  50  0  1  1  1  1  -1   0\n"
                           b"        10000   9.447840918743423E+03 "
                           b"6.36184172E-02 /A\n"
                           b"        20000   9.246852060423453E+03 "
                           b"1.63997319E-02 /B\n"
                           b"        30000   4.671341872421471E+03 "
                           b"5.27198401E-02 /C")
            var_file.seek(0)
            fit_file.write(b'I am a fit file')
            fit_file.seek(0)
            lin_file.write(b'I am a lin file')
            lin_file.seek(0)
            qpart_file.write(b"#form : interpolation\n"
                             b"200.000   331777.6674")
            qpart_file.seek(0)
            payload = {
                'species': species.id,
                'molecule_tag': 1,
                'hyperfine': True,
                'degree_of_freedom': 3,
                'category': 'asymmetric top',
                'partition_function': json.dumps({'300.000': '331777.6674'}),
                'mu_a': 0.5,
                'mu_b': 0,
                'mu_c': 0,
                'a_const': 0.123,
                'b_const': 0,
                'c_const': 0,
                'linelist': linelist.id,
                'data_date': '2020-01-01',
                'data_contributor': 'Test Contributor',
                'int_file': int_file,
                'var_file': var_file,
                'fit_file': fit_file,
                'lin_file': lin_file,
                'qpart_file': qpart_file,
                'notes': 'Test Species Metadata Put',
            }
            res = self.client.put(url, payload, format='multipart')
        self.assertEqual(res.status_code, status.HTTP_400_BAD_REQUEST)

    def test_delete_meta(self):
        """Test deleting a species metadata."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        url = reverse('data:speciesmetadata-detail',
                      args=[meta.id]) + '?delete_reason=Test delete reason'

        res = self.client.delete(url)
        self.assertEqual(res.status_code, status.HTTP_204_NO_CONTENT)
        self.assertFalse(SpeciesMetadata.objects.filter(id=meta.id).exists())
        self.assertEqual(SpeciesMetadata.history.filter(id=meta.id).first(
        ).history_change_reason, 'Test delete reason')
        self.assertEqual(SpeciesMetadata.history.filter(id=meta.id).first(
        ).history_user_id, self.user.id)
        history_count = SpeciesMetadata.history.filter(
            id=meta.id).count()
        self.assertEqual(history_count, 2)

    def test_delete_meta_without_delete_reason_fails(self):
        """Test deleting a species metadata without delete reason fails."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        url = reverse('data:speciesmetadata-detail', args=[meta.id])

        res = self.client.delete(url)
        self.assertEqual(res.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertTrue(SpeciesMetadata.objects.filter(id=meta.id).exists())

    def test_delete_referenced_meta_fails(self):
        """Test deleting a species metadata that is referenced fails."""
        species = create_species()
        linelist = create_linelist()
        meta = create_meta(species.id, linelist.id)
        create_line(meta.id)
        url = reverse('data:speciesmetadata-detail',
                      args=[meta.id]) + '?delete_reason=Test delete reason'

        res = self.client.delete(url)
        self.assertEqual(res.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertTrue(SpeciesMetadata.objects.filter(id=meta.id).exists())
