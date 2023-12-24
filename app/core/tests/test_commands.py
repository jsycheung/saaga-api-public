"""
Test custom Django management commands
"""

# SimpleTestCase is used because we do not need migrations to be applied.
# We are just simulating behavior of the database.
# We use SimpleTestCase so that it does not create any database setup.
from django.test import SimpleTestCase

# patch allows mocking of behavior of the database
# to simulate when the database is returning a response or not.
from unittest.mock import patch

from psycopg2 import OperationalError as Psycopg2Error
from django.db.utils import OperationalError
from django.core.management import call_command

"""
Patch decorator mocks the command 'check', or the behavior of the database.
'check' is provided by the BaseCommand class.
'check' is a method that checks the status of the database.
Here we simulate the response of the 'check' method.
Patch adds a new argument to each of the calls as patched_check.
"""


@patch('core.management.commands.wait_for_db.Command.check')
class CommandTests(SimpleTestCase):
    """Test commands."""

    def test_wait_for_db_ready(self, patched_check):
        """Test waiting for database if database is ready."""
        # When check is called inside the test case, it returns the true value.
        patched_check.return_value = True
        # Check if the command can be called.
        call_command('wait_for_db')
        # Check if the check method is called with correct databases parameter.
        patched_check.assert_called_once_with(databases=['default'])

    @patch('time.sleep')
    # Replace the sleep function and just return a None value.
    # We override the behavior of sleep so it does not wait during testing.
    def test_wait_for_db_delay(self, patched_sleep, patched_check):
        """
        Test waiting for database when getting OperationalError.
        Make the check method raise exceptions -
        Psycopg2Error for the first two times,
        OperationalError for the next three times,
        and return true value at last.
        """
        patched_check.side_effect = [Psycopg2Error] * 2 + \
            [OperationalError] * 3 + [True]

        # When wait_for_db command is called, the check method raise exceptions
        # or return values according to patched_check.side_effect.
        call_command('wait_for_db')

        # Make sure the check method is called six times.
        self.assertEqual(patched_check.call_count, 6)

        # Check if the check method is called with correct databases parameter.
        patched_check.assert_called_with(databases=['default'])
