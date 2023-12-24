"""
URL mappings for the data app.
"""
from django.urls import (
    path,
    include,
)

from rest_framework.routers import DefaultRouter

from data import views

router = DefaultRouter()
router.register('species', views.SpeciesViewSet)
router.register('linelist', views.LinelistViewSet)
router.register('species-metadata', views.SpeciesMetadataViewSet)
router.register('reference', views.ReferenceViewSet)
router.register('meta-reference', views.MetaReferenceViewSet)
router.register('line', views.LineViewSet)

app_name = 'data'

urlpatterns = [
    path('', include(router.urls))
]
