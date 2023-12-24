"""
Test for models
"""
from django.test import TestCase
from django.contrib.auth import get_user_model
from core import models
import datetime
from django.core.files.uploadedfile import SimpleUploadedFile
from decimal import Decimal
from rdkit import Chem
from rdkit.Chem import Descriptors
import selfies as sf
from django.db.utils import IntegrityError


def create_user(email="test@example.com", password="testpw123",
                name="test person", organization="test org"):
    """Create and return a new user."""
    return get_user_model().objects.create_user(
        email=email,
        password=password,
        name=name,
        organization=organization
    )


def create_species(standard_inchi="test standard inchi",
                   name=["test species 1", "test species 2"],
                   iupac_name="test iupac name",
                   name_formula="test name formula",
                   name_html="test name html",
                   smiles="C",
                   standard_inchi_key="test standard inchi key 1",
                   notes="test"):
    """Create and return a new species."""
    return models.Species.objects.create(
        name=name,
        iupac_name=iupac_name,
        name_formula=name_formula,
        name_html=name_html,
        molecular_mass=Descriptors.ExactMolWt(Chem.MolFromSmiles(smiles)),
        smiles=smiles,
        standard_inchi=standard_inchi,
        standard_inchi_key=standard_inchi_key,
        selfies=sf.encoder(smiles),
        mol_obj=Chem.MolFromSmiles(smiles),
        notes=notes
    )


def create_linelist(linelist_name="test linelist"):
    """Create and return a new linelist."""
    return models.Linelist.objects.create(
        linelist_name=linelist_name
    )


def create_metadata(species, linelist,
                    molecule_tag=1,
                    hyperfine=True,
                    degree_of_freedom=3,
                    category="test category",
                    partition_function={"300": "1.0"},
                    mu_a=Decimal("1.2"), mu_b=Decimal("1.2"),
                    mu_c=Decimal("1.2"),
                    a_const=Decimal("1.2"),
                    b_const=Decimal("1.2"),
                    c_const=Decimal("1.2"),
                    data_date=datetime.date.today(),
                    data_contributor="test data contributor",
                    int_file="test_int_file",
                    var_file="test_var_file",
                    fit_file="test_fit_file",
                    lin_file="test_lin_file",
                    qpart_file="test_qpart_file",
                    notes="test"):
    return models.SpeciesMetadata.objects.create(
        species=species,
        molecule_tag=molecule_tag,
        hyperfine=hyperfine,
        degree_of_freedom=degree_of_freedom,
        category=category,
        partition_function=partition_function,
        mu_a=mu_a,
        mu_b=mu_b,
        mu_c=mu_c,
        a_const=a_const,
        b_const=b_const,
        c_const=c_const,
        linelist=linelist,
        data_date=data_date,
        data_contributor=data_contributor,
        int_file=int_file,
        var_file=var_file,
        fit_file=fit_file,
        lin_file=lin_file,
        qpart_file=qpart_file,
        notes=notes
    )


def create_reference(doi="test doi",
                     ref_url="test url",
                     bibtex="bibtex url",
                     notes=""):
    return models.Reference.objects.create(
        doi=doi,
        ref_url=ref_url,
        bibtex=bibtex,
        notes=notes
    )


class ModelTests(TestCase):
    """Test models."""

    def test_create_user_successful(self):
        """Test creating a user with an email is successful."""
        email = "test@example.com"
        password = "testpw123"
        name = "test person"
        organization = "test org"
        user = get_user_model().objects.create_user(
            email=email,
            password=password,
            name=name,
            organization=organization
        )
        self.assertEqual(user.email, email)
        self.assertTrue(user.check_password(password))
        self.assertEqual(user.name, name)
        self.assertEqual(user.organization, organization)
        self.assertTrue(user.is_active)
        self.assertTrue(user.is_staff)
        self.assertTrue(user.is_superuser)

    def test_new_user_email_normalized(self):
        """Test email is normalized for new users."""
        sample_emails = [
            ['test1@EXAMPLE.com', 'test1@example.com'],
            ['Test2@Example.com', 'Test2@example.com'],
            ['TEST3@EXAMPLE.COM', 'TEST3@example.com'],
            ['test4@example.com', 'test4@example.com']
        ]
        for email, expected in sample_emails:
            user = get_user_model().objects.create_user(
                email, 'pw123', 'name123', 'org123')
            self.assertEqual(user.email, expected)

    def test_new_user_no_required_fields_raises_error(self):
        """Test that creating a user without email, name,
        or organization raises an error."""
        with self.assertRaises(ValueError):
            get_user_model().objects.create_user(
                '', 'pw123', 'name123', 'org123')
        with self.assertRaises(ValueError):
            get_user_model().objects.create_user(
                'test@example.com', '', 'name123', 'org123')
        with self.assertRaises(ValueError):
            get_user_model().objects.create_user(
                'test@example.com', 'pw123', '', 'org123')
        with self.assertRaises(ValueError):
            get_user_model().objects.create_user(
                'test@example.com', 'pw123', 'name123', '')

    def test_create_superuser(self):
        """Test creating a superuser."""
        email = "test@example.com"
        password = "testpw123"
        name = "test person"
        organization = "test org"
        user = get_user_model().objects.create_superuser(
            email=email,
            password=password,
            name=name,
            organization=organization
        )
        self.assertTrue(user.is_active)
        self.assertTrue(user.is_staff)
        self.assertTrue(user.is_superuser)

    def test_create_user_duplicate(self):
        """Test creating a user with duplicate email raises an error."""
        create_user(email="testuser@example.com")
        with self.assertRaises(IntegrityError):
            create_user(email="testuser@example.com")

    def test_create_linelist_lower_case(self):
        """Test creating a linelist with lower case name is successful."""
        linelist = models.Linelist.objects.create(
            linelist_name="test linelist"
        )
        self.assertEqual(str(linelist), linelist.linelist_name)
        self.assertEqual(linelist.linelist_name, "test linelist")

    def test_create_linelist_upper_case(self):
        """Test creating a linelist with upper case name is successful."""
        linelist = models.Linelist.objects.create(
            linelist_name="TEST LINELIST"
        )
        self.assertEqual(str(linelist), linelist.linelist_name)
        self.assertEqual(linelist.linelist_name, "test linelist")

    def test_create_linelist_duplicate(self):
        """Test creating a linelist with duplicate name raises an error."""
        create_linelist(linelist_name="test linelist duplicate")
        with self.assertRaises(IntegrityError):
            create_linelist(linelist_name="test linelist duplicate")

    def test_create_reference(self):
        """Test creating a reference is successful."""
        reference = models.Reference.objects.create(
            doi="10.3847/1538-4357/acc584",
            ref_url="https://iopscience.iop.org/article/ \
            10.3847/1538-4357/acc584/pdf",
            bibtex=SimpleUploadedFile(
                "test.bib",
                b"@ARTICLE{2023ApJ...948..133C,}"
            ),
            notes=""
        )
        self.assertEqual(str(reference), reference.ref_url)

    def test_create_reference_duplicate(self):
        """Test creating a reference with duplicate ref_url raises an error."""
        create_reference(ref_url='123')
        with self.assertRaises(IntegrityError):
            create_reference(ref_url='123')

    def test_create_species(self):
        """Test creating a species is successful."""
        species = models.Species.objects.create(
            name="test species",
            iupac_name="test iupac name",
            name_formula="test name formula",
            name_html="test name html",
            smiles="CC",
            molecular_mass=Descriptors.ExactMolWt(Chem.MolFromSmiles("CC")),
            standard_inchi="test standard inchi",
            standard_inchi_key="test standard inchi key",
            selfies=sf.encoder("CC"),
            mol_obj=Chem.MolFromSmiles("CC"),
            notes="test"
        )
        self.assertEqual(str(species), species.iupac_name)

    def test_create_species_duplicate(self):
        """Test creating a species with
        duplicate standard_inchi raises error."""
        create_species(standard_inchi="test standard inchi duplicate")
        with self.assertRaises(IntegrityError):
            create_species(standard_inchi="test standard inchi duplicate")

    def test_create_speciesmetadata(self):
        """Test creating a species metadata is successful."""
        species = create_species(
            standard_inchi="test standard inchi metadata",)
        linelist = create_linelist(linelist_name="test linelist metadata")
        speciesmetadata = models.SpeciesMetadata.objects.create(
            species=species,
            molecule_tag=1,
            hyperfine=True,
            degree_of_freedom=3,
            category="test category",
            partition_function={"300": "1.0"},
            mu_a=Decimal("1.2"),
            mu_b=Decimal("1.2"),
            mu_c=Decimal("1.2"),
            a_const=Decimal("1.2"),
            b_const=Decimal("1.2"),
            c_const=Decimal("1.2"),
            linelist=linelist,
            data_date=datetime.date.today(),
            data_contributor="test data contributor",
            int_file="test_int_file",
            var_file="test_var_file",
            fit_file="test_fit_file",
            lin_file="test_lin_file",
            qpart_file="test_qpart_file",
            notes="test"
        )
        self.assertEqual(str(speciesmetadata),
                         "species metadata of " +
                         speciesmetadata.species.iupac_name)

    def test_create_metareference_dipole(self):
        """Test creating a metadata reference
        for dipole moment is successful."""
        linelist = create_linelist(linelist_name="test linelist dipole moment")
        species = create_species(
            standard_inchi="test standard inchi dipole moment")
        metadata = create_metadata(linelist=linelist, species=species)
        ref = create_reference()
        metaref = models.MetaReference.objects.create(
            meta=metadata,
            ref=ref,
            dipole_moment=True,
            spectrum=False,
            notes="test notes"
        )
        self.assertEqual(str(
            metaref), "metadata reference for dipole moment of " +
            metaref.meta.species.iupac_name)

    def test_create_metareference_spectrum(self):
        """Test creating a metadata reference for spectrum is successful."""
        linelist = create_linelist(linelist_name="test linelist spectrum")
        species = create_species(standard_inchi="test standard inchi spectrum")
        metadata = create_metadata(linelist=linelist, species=species)
        ref = create_reference()
        metaref = models.MetaReference.objects.create(
            meta=metadata,
            ref=ref,
            dipole_moment=False,
            spectrum=True,
            notes="test notes"
        )
        self.assertEqual(str(
            metaref), "metadata reference for spectrum of " +
            metaref.meta.species.iupac_name)

    def test_create_metareference_dipole_spectrum(self):
        """Test creating a metadata reference for
        dipole moment and spectrum is successful."""
        linelist = create_linelist(
            linelist_name="test linelist dipole spectrum")
        species = create_species(
            standard_inchi="test standard inchi dipole spectrum")
        metadata = create_metadata(linelist=linelist, species=species)
        ref = create_reference()
        metaref = models.MetaReference.objects.create(
            meta=metadata,
            ref=ref,
            dipole_moment=True,
            spectrum=True,
            notes="test notes"
        )
        self.assertEqual(str(
            metaref), "metadata reference for dipole moment and spectrum of " +
            metaref.meta.species.iupac_name)

    def test_create_line(self):
        """Test creating a line is successful."""
        linelist = create_linelist(linelist_name="test linelist line")
        species = create_species(
            standard_inchi="test standard inchi line")
        metadata = create_metadata(linelist=linelist, species=species)
        line = models.Line.objects.create(
            meta=metadata,
            measured=True,
            frequency=1234.56,
            uncertainty=0.12,
            intensity=5678.1234,
            s_ij=123.456,
            s_ij_mu2=678567.167,
            a_ij=67867.1267,
            lower_state_energy=123456778.145,
            upper_state_energy=1267890.1345,
            lower_state_degeneracy=2,
            upper_state_degeneracy=3,
            lower_state_qn={"test": "qn"},
            upper_state_qn={"test": "qn"},
            rovibrational=True,
            vib_qn="test vib qn",
            pickett_qn_code=123,
            pickett_lower_state_qn="010101",
            pickett_upper_state_qn="010102",
            notes="test notes"
        )
        self.assertEqual(str(line), "line of "+line.meta.species.iupac_name)
