"""
Tests for the Django admin modification.
"""
from django.test import TestCase, Client
from django.contrib.auth import get_user_model
from django.urls import reverse


class AdminSiteTests(TestCase):
    """Tests for Django admin."""

    def setUp(self):
        """"Create user and client."""
        self.client = Client()
        self.admin_user = get_user_model().objects.create_superuser(
            email="admin@example.com",
            password="adminpw123",
            name="test person admin",
            organization="test org admin"
        )
        self.client.force_login(self.admin_user)
        self.user = get_user_model().objects.create_user(
            email="test@example.com",
            password="testpw123",
            name="test person",
            organization="test org"
        )

    def test_users_list(self):
        """Test that users are listed on page."""
        url = reverse("admin:core_user_changelist")
        res = self.client.get(url)
        self.assertContains(res, self.admin_user.email)
        self.assertContains(res, self.admin_user.name)
        self.assertContains(res, self.admin_user.organization)
        self.assertContains(res, self.user.email)
        self.assertContains(res, self.user.name)
        self.assertContains(res, self.user.organization)

    def test_edit_user_page(self):
        """Test the edit user page works."""
        url = reverse('admin:core_user_change', args=[self.user.id])
        res = self.client.get(url)
        self.assertEqual(res.status_code, 200)

    def test_create_user_page(self):
        """Test the create user page works."""
        url = reverse('admin:core_user_add')
        res = self.client.get(url)

        self.assertEqual(res.status_code, 200)
