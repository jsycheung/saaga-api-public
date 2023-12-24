"""
Serializers for data APIs.
"""
from rest_framework import serializers

from core.models import (Species, Linelist, SpeciesMetadata,
                         Reference, MetaReference, Line)


class LinelistSerializer(serializers.ModelSerializer):
    """Serialzer for linelists."""
    class Meta:
        model = Linelist
        fields = ['id', 'linelist_name']
        read_only_fields = ['id']


class LinelistChangeSerializer(LinelistSerializer):
    """Serializer for put and patch linelists."""
    _change_reason = serializers.CharField(
        max_length=255, write_only=True, required=True)

    class Meta(LinelistSerializer.Meta):
        fields = LinelistSerializer.Meta.fields + ['_change_reason']


class ReferenceSerializer(serializers.ModelSerializer):
    """Serializer for references."""
    class Meta:
        model = Reference
        fields = ['id', 'doi', 'ref_url', 'bibtex',
                  'notes']
        read_only_fields = ['id']


class ReferenceChangeSerializer(ReferenceSerializer):
    """Serializer for put and patch references."""
    _change_reason = serializers.CharField(
        max_length=255, write_only=True, required=True)

    class Meta(ReferenceSerializer.Meta):
        fields = ReferenceSerializer.Meta.fields + ['_change_reason']


class SpeciesSerializer(serializers.ModelSerializer):
    """Serializer for species."""
    molecular_mass = serializers.DecimalField(
        max_digits=None, decimal_places=None, read_only=True)

    class Meta:
        model = Species
        fields = ['id', 'name', 'iupac_name', 'name_formula',
                  'name_html', 'molecular_mass', 'smiles',
                  'standard_inchi', 'standard_inchi_key', 'selfies', 'notes']
        read_only_fields = ['id', 'molecular_mass', 'selfies']


class SpeciesChangeSerializer(SpeciesSerializer):
    """Serializer for put and patch species."""
    _change_reason = serializers.CharField(
        max_length=255, write_only=True, required=True)

    class Meta(SpeciesSerializer.Meta):
        fields = SpeciesSerializer.Meta.fields + ['_change_reason', 'mol_obj']
        read_only_fields = ['id']


class SpeciesMetadataSerializer(serializers.ModelSerializer):
    """Serializer for species metadata."""
    mu_a = serializers.DecimalField(
        max_digits=None, decimal_places=None, required=False, allow_null=True)
    mu_b = serializers.DecimalField(
        max_digits=None, decimal_places=None, required=False, allow_null=True)
    mu_c = serializers.DecimalField(
        max_digits=None, decimal_places=None, required=False, allow_null=True)
    a_const = serializers.DecimalField(
        max_digits=None, decimal_places=None, required=False, allow_null=True)
    b_const = serializers.DecimalField(
        max_digits=None, decimal_places=None, required=False, allow_null=True)
    c_const = serializers.DecimalField(
        max_digits=None, decimal_places=None, required=False, allow_null=True)

    class Meta:
        model = SpeciesMetadata
        fields = ['id', 'species', 'molecule_tag', 'hyperfine',
                  'degree_of_freedom', 'category', 'partition_function',
                  'mu_a', 'mu_b', 'mu_c', 'a_const', 'b_const', 'c_const',
                  'linelist', 'data_date', 'data_contributor', 'qpart_file',
                  'int_file', 'var_file', 'fit_file', 'lin_file', 'notes']
        read_only_fields = ['id', 'partition_function']


class SpeciesMetadataChangeSerializer(SpeciesMetadataSerializer):
    """Serializer for put and patch species metadata."""
    _change_reason = serializers.CharField(
        max_length=255, write_only=True, required=True)

    class Meta(SpeciesMetadataSerializer.Meta):
        fields = SpeciesMetadataSerializer.Meta.fields + ['_change_reason']
        read_only_fields = ['id']


class MetaReferenceSerializer(serializers.ModelSerializer):
    """Serializer for metadata references."""
    class Meta:
        model = MetaReference
        fields = ['id', 'meta', 'ref', 'dipole_moment',
                  'spectrum', 'notes']
        read_only_fields = ['id']


class MetaReferenceChangeSerializer(MetaReferenceSerializer):
    """Serializer for put and patch metadata references."""
    _change_reason = serializers.CharField(
        max_length=255, write_only=True, required=True)

    class Meta(MetaReferenceSerializer.Meta):
        fields = MetaReferenceSerializer.Meta.fields + ['_change_reason']


class LineSerializer(serializers.ModelSerializer):
    """Serializer for creating lines with POST request."""
    frequency = serializers.DecimalField(
        max_digits=None, decimal_places=None, read_only=True)
    uncertainty = serializers.DecimalField(
        max_digits=None, decimal_places=None, read_only=True)
    intensity = serializers.DecimalField(
        max_digits=None, decimal_places=None, read_only=True)
    s_ij = serializers.DecimalField(
        max_digits=None, decimal_places=None, read_only=True)
    s_ij_mu2 = serializers.DecimalField(
        max_digits=None, decimal_places=None, read_only=True)
    a_ij = serializers.DecimalField(
        max_digits=None, decimal_places=None, read_only=True)
    lower_state_energy = serializers.DecimalField(
        max_digits=None, decimal_places=None, read_only=True)
    upper_state_energy = serializers.DecimalField(
        max_digits=None, decimal_places=None, read_only=True)
    cat_file = serializers.FileField(write_only=True)
    qn_label_str = serializers.CharField(write_only=True)
    contains_rovibrational = serializers.BooleanField(write_only=True)

    class Meta:
        model = Line
        fields = ['id', 'meta', 'measured', 'cat_file', 'qn_label_str',
                  'frequency', 'uncertainty', 'intensity',
                  's_ij', 's_ij_mu2', 'a_ij', 'upper_state_energy',
                  'lower_state_energy', 'upper_state_degeneracy',
                  'lower_state_degeneracy', 'upper_state_qn',
                  'lower_state_qn', 'contains_rovibrational',
                  'rovibrational', 'vib_qn', 'pickett_qn_code',
                  'pickett_upper_state_qn', 'pickett_lower_state_qn',
                  'notes']
        read_only_fields = ['id', 'measured', 'frequency', 'uncertainty',
                            'intensity', 's_ij', 's_ij_mu2', 'a_ij',
                            'lower_state_energy', 'upper_state_energy',
                            'lower_state_degeneracy',
                            'upper_state_degeneracy', 'upper_state_qn',
                            'lower_state_qn', 'rovibrational',
                            'pickett_qn_code', 'pickett_upper_state_qn',
                            'pickett_lower_state_qn']


class LineSerializerList(serializers.ModelSerializer):
    """Serializer for creating lines in the backend
    after receiving POST request."""
    frequency = serializers.DecimalField(
        max_digits=None, decimal_places=None)
    uncertainty = serializers.DecimalField(
        max_digits=None, decimal_places=None)
    intensity = serializers.DecimalField(
        max_digits=None, decimal_places=None)
    s_ij = serializers.DecimalField(
        max_digits=None, decimal_places=None, required=False, allow_null=True)
    s_ij_mu2 = serializers.DecimalField(
        max_digits=None, decimal_places=None)
    a_ij = serializers.DecimalField(
        max_digits=None, decimal_places=None)
    lower_state_energy = serializers.DecimalField(
        max_digits=None, decimal_places=None)
    upper_state_energy = serializers.DecimalField(
        max_digits=None, decimal_places=None)

    class Meta:
        model = Line
        fields = ['id', 'meta', 'measured', 'frequency', 'uncertainty',
                  'intensity', 's_ij', 's_ij_mu2', 'a_ij',
                  'upper_state_energy', 'lower_state_energy',
                  'upper_state_degeneracy', 'lower_state_degeneracy',
                  'upper_state_qn', 'lower_state_qn', 'rovibrational',
                  'vib_qn', 'pickett_qn_code', 'pickett_upper_state_qn',
                  'pickett_lower_state_qn', 'notes']
        read_only_fields = ['id']


class LineChangeSerializerList(LineSerializerList):
    """Serializer for put and patch lines."""
    _change_reason = serializers.CharField(
        max_length=255, write_only=True, required=True)

    class Meta(LineSerializerList.Meta):
        fields = LineSerializerList.Meta.fields + ['_change_reason']


class QuerySerializer(serializers.ModelSerializer):
    """Serilizer for querying lines falling
    within specified frequency range."""
    frequency = serializers.DecimalField(max_digits=None, decimal_places=None)
    uncertainty = serializers.DecimalField(
        max_digits=None, decimal_places=None)
    intensity = serializers.DecimalField(max_digits=None, decimal_places=None)
    s_ij = serializers.DecimalField(max_digits=None, decimal_places=None)
    s_ij_mu2 = serializers.DecimalField(max_digits=None, decimal_places=None)
    a_ij = serializers.DecimalField(max_digits=None, decimal_places=None)
    lower_state_energy = serializers.DecimalField(
        max_digits=None, decimal_places=None)
    upper_state_energy = serializers.DecimalField(
        max_digits=None, decimal_places=None)

    class Meta:
        model = Line
        fields = ['frequency', 'measured', 'uncertainty', 'intensity',
                  'lower_state_qn', 'upper_state_qn',
                  'lower_state_energy', 'upper_state_energy',
                  's_ij', 's_ij_mu2', 'a_ij',
                  'rovibrational']

    def to_representation(self, instance):
        representation = super().to_representation(instance)
        representation['name_formula'] = instance.meta.species.name_formula
        representation['iupac_name'] = instance.meta.species.iupac_name
        representation['name'] = instance.meta.species.name
        representation['molecule_tag'] = instance.meta.molecule_tag
        representation['hyperfine'] = instance.meta.hyperfine
        representation['linelist'] = instance.meta.linelist.linelist_name
        representation['meta_id'] = instance.meta.id
        representation['smiles'] = instance.meta.species.smiles
        representation['selfies'] = instance.meta.species.selfies
        return representation
