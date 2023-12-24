"""
Django admin customization.
"""

from django.contrib import admin
from django.contrib.auth.admin import UserAdmin as BaseUserAdmin
from django.utils.translation import gettext_lazy as _

from core import models


class UserAdmin(BaseUserAdmin):
    """Define the admin pages for users."""
    ordering = ['id']
    list_display = ['email', 'name', 'organization']
    fieldsets = (
        (None, {'fields': ('email', 'password', 'name', 'organization')}),
        (
            _('Permissions'),
            {
                'fields': (
                    'is_active',
                    'is_staff',
                    'is_superuser'
                )
            }
        ),
        (_('Important dates'), {'fields': ('last_login',)})
    )
    readonly_fields = ['last_login']
    add_fieldsets = (
        (None, {
            'classes': ('wide',),
            'fields': (
                'email',
                'password1',
                'password2',
                'name',
                'organization',
                'is_active',
                'is_staff',
                'is_superuser'
            )
        }),
    )


class SpeciesAdmin(admin.ModelAdmin):
    """Define the admin pages for users."""
    ordering = ['id']
    list_display = ['iupac_name', 'name_formula', 'smiles']
    fieldsets = (
        (_('Names'), {
         'fields': ('name', 'iupac_name', 'name_formula', 'name_html')}),
        (
            _('Identifiers'),
            {
                'fields': (
                    'smiles',
                    'standard_inchi',
                    'standard_inchi_key',
                    'selfies',
                    'display_mol',
                )
            }
        ),
        (_('Records'), {'fields': ('notes',)})
    )
    readonly_fields = ['selfies', 'display_mol']
    add_fieldsets = (
        (None, {
            'classes': ('wide',),
            'fields': (
                'name',
                'iupac_name',
                'name_formula',
                'name_html',
                'canonical_smiles',
                'standard_inchi',
                'standard_inchi_key',
                'notes'
            )
        }),
    )


admin.site.register(models.User, UserAdmin)
admin.site.register(models.Line)
admin.site.register(models.SpeciesMetadata)
admin.site.register(models.MetaReference)
admin.site.register(models.Reference)
admin.site.register(models.Linelist)
admin.site.register(models.Species, SpeciesAdmin)
