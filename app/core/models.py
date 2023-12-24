"""
Database models.
"""
import os
import uuid
from django_rdkit import models
from django.contrib.auth.models import (
    AbstractBaseUser,
    BaseUserManager,
    PermissionsMixin
)
from django.utils.functional import cached_property
from django.utils.html import format_html
from rdkit.Chem import Draw
import base64
from django.core.validators import FileExtensionValidator
from simple_history.models import HistoricalRecords
from simple_history import register
from django.contrib.postgres.indexes import GistIndex


class ArbitraryDecimalField(models.DecimalField):
    """
    Custom DecimalField that allows for arbitrary precision
    (max_digits) and scale (decimal_places).
    """

    def _check_decimal_places(self):
        return []

    def _check_max_digits(self):
        return []

    def _check_decimal_places_and_max_digits(self, **kwargs):
        return []

    def db_type(self, connection):
        # pg or bust
        assert connection.settings_dict["ENGINE"] == \
            "django.db.backends.postgresql"
        return "numeric"


def sp_file_path(instance, filename):
    """Generate file path for SPFIT/SPCAT
    (.int, .var, .lin, .fit, .qpart) files."""
    ext = os.path.splitext(filename)[1]
    filename = f'{uuid.uuid4()}{ext}'
    return os.path.join('uploads', 'sp', filename)


def bib_file_path(instance, filename):
    """Generage file path for .bib files"""
    ext = os.path.splitext(filename)[1]
    filename = f'{uuid.uuid4()}{ext}'
    return os.path.join('uploads', 'bib', filename)


class UserManager(BaseUserManager):
    """Manager for users."""

    def create_user(self, email, password, name, organization):
        """Create, save and return a new user."""
        if not email or not password or not name or not organization:
            raise ValueError('User must have an email, \
                              password, name, and organization')
        user = self.model(email=self.normalize_email(email),
                          name=name, organization=organization)
        user.set_password(password)
        user.save(using=self._db)

        return user

    def create_superuser(self, email, password, name, organization):
        """Create and return a new superuser."""
        user = self.create_user(email, password, name, organization)
        user.is_staff = True
        user.is_superuser = True
        user.save(using=self._db)

        return user


class User(AbstractBaseUser, PermissionsMixin):
    """Users in the system."""
    email = models.EmailField(max_length=255, unique=True)
    name = models.CharField(max_length=255)
    organization = models.CharField(max_length=255)
    is_active = models.BooleanField(default=True)
    is_staff = models.BooleanField(default=True)
    is_superuser = models.BooleanField(default=True)

    REQUIRED_FIELDS = ['name', 'organization']
    objects = UserManager()

    USERNAME_FIELD = 'email'


register(User)


class Linelist(models.Model):
    """Linelist object."""
    linelist_name = models.CharField(max_length=255, unique=True)
    history = HistoricalRecords()

    def save(self, *args, **kwargs):
        self.linelist_name = self.linelist_name.lower()
        return super(Linelist, self).save(*args, **kwargs)

    def __str__(self):
        return self.linelist_name


class Reference(models.Model):
    """Reference object."""
    doi = models.CharField(max_length=255, blank=True)
    ref_url = models.CharField(max_length=255, unique=True)
    bibtex = models.FileField(upload_to=bib_file_path,
                              validators=[FileExtensionValidator(
                                  allowed_extensions=["bib"])])
    notes = models.TextField(blank=True)
    history = HistoricalRecords()

    def __str__(self):
        return self.ref_url


class Species(models.Model):
    """Species object."""
    class Meta:
        verbose_name_plural = 'Species'
        indexes = [
            GistIndex(fields=['mol_obj']),
        ]
    name = models.JSONField()
    iupac_name = models.CharField(max_length=255, unique=True)
    name_formula = models.CharField(max_length=255)
    name_html = models.CharField(max_length=255)
    molecular_mass = ArbitraryDecimalField()
    smiles = models.CharField(max_length=255)
    standard_inchi = models.CharField(max_length=255)
    standard_inchi_key = models.CharField(max_length=255)
    selfies = models.CharField(max_length=255)
    mol_obj = models.MolField()
    notes = models.TextField(blank=True)
    history = HistoricalRecords()

    def __str__(self):
        return self.iupac_name

    @cached_property
    def display_mol(self):
        """function for displaying the rdkit mol object in the
        form of image in django admin."""
        if self.mol_obj:
            dm = Draw.PrepareMolForDrawing(self.mol_obj)
            d2d = Draw.MolDraw2DCairo(400, 400)
            d2d.DrawMolecule(dm)
            d2d.FinishDrawing()
            text = d2d.GetDrawingText()
            imtext = base64.b64encode(text).decode('utf8')
            html = '<img src="data:image/png;base64, {img}" alt="rdkit image">'
            return format_html(html, img=imtext)
        return format_html('<strong>There is no image for this entry.<strong>')
    display_mol.short_description = 'Display rdkit image'


class SpeciesMetadata(models.Model):
    """Species metadata object."""
    class Meta:
        verbose_name_plural = 'Species metadata'
    species = models.ForeignKey(
        'Species',
        on_delete=models.PROTECT
    )
    molecule_tag = models.IntegerField(blank=True, null=True)
    hyperfine = models.BooleanField()
    degree_of_freedom = models.IntegerField()
    category = models.CharField(max_length=255)
    partition_function = models.JSONField()
    mu_a = ArbitraryDecimalField(null=True)
    mu_b = ArbitraryDecimalField(null=True)
    mu_c = ArbitraryDecimalField(null=True)
    a_const = ArbitraryDecimalField(null=True)
    b_const = ArbitraryDecimalField(null=True)
    c_const = ArbitraryDecimalField(null=True)
    linelist = models.ForeignKey(
        'Linelist',
        on_delete=models.PROTECT
    )
    data_date = models.DateField()
    data_contributor = models.CharField(max_length=255)
    int_file = models.FileField(null=True, upload_to=sp_file_path,
                                validators=[FileExtensionValidator(
                                    allowed_extensions=["int"])])
    var_file = models.FileField(null=True, upload_to=sp_file_path,
                                validators=[FileExtensionValidator(
                                    allowed_extensions=["var"])])
    fit_file = models.FileField(null=True, upload_to=sp_file_path,
                                validators=[FileExtensionValidator(
                                    allowed_extensions=["fit"])])
    lin_file = models.FileField(null=True, upload_to=sp_file_path,
                                validators=[FileExtensionValidator(
                                    allowed_extensions=["lin"])])
    qpart_file = models.FileField(upload_to=sp_file_path, validators=[
                                  FileExtensionValidator(
                                      allowed_extensions=["qpart"])])
    notes = models.TextField(blank=True)
    history = HistoricalRecords()

    def __str__(self):
        return "species metadata of "+self.species.iupac_name


class MetaReference(models.Model):
    """Metadata reference object relating species metadata with references"""
    meta = models.ForeignKey(
        'SpeciesMetadata',
        on_delete=models.PROTECT
    )
    ref = models.ForeignKey(
        'Reference',
        on_delete=models.PROTECT
    )
    dipole_moment = models.BooleanField()
    spectrum = models.BooleanField()
    notes = models.TextField(blank=True)
    history = HistoricalRecords()

    def __str__(self):
        if self.dipole_moment and self.spectrum:
            return "metadata reference for dipole moment "\
                "and spectrum of " + self.meta.species.iupac_name
        elif self.dipole_moment:
            return "metadata reference for dipole moment of " + \
                self.meta.species.iupac_name
        elif self.spectrum:
            return "metadata reference for spectrum of " + \
                self.meta.species.iupac_name
        else:
            return "metadata reference for "+self.meta.species.iupac_name


class Line(models.Model):
    """Line object."""
    meta = models.ForeignKey(
        'SpeciesMetadata',
        on_delete=models.PROTECT
    )
    measured = models.BooleanField()
    frequency = ArbitraryDecimalField()
    uncertainty = ArbitraryDecimalField()
    intensity = ArbitraryDecimalField()
    s_ij = ArbitraryDecimalField(null=True)
    s_ij_mu2 = ArbitraryDecimalField()
    a_ij = ArbitraryDecimalField()
    lower_state_energy = ArbitraryDecimalField()
    upper_state_energy = ArbitraryDecimalField()
    lower_state_degeneracy = models.IntegerField()
    upper_state_degeneracy = models.IntegerField()
    lower_state_qn = models.JSONField()
    upper_state_qn = models.JSONField()
    rovibrational = models.BooleanField()
    vib_qn = models.CharField(max_length=255, blank=True)
    pickett_qn_code = models.IntegerField()
    pickett_lower_state_qn = models.CharField(max_length=255)
    pickett_upper_state_qn = models.CharField(max_length=255)
    notes = models.TextField(blank=True)
    history = HistoricalRecords()

    class Meta:
        indexes = [
            models.Index(fields=['frequency'])
        ]

    def __str__(self):
        return "line of "+self.meta.species.iupac_name
