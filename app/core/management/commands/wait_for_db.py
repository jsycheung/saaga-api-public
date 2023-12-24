"""
Django command to wait for the database to be available.
"""
import time
from psycopg2 import OperationalError as Psycopg2OpError
from django.db.utils import OperationalError
from django.core.management.base import BaseCommand


class Command(BaseCommand):
    """Django command to wait for database"""

    def handle(self, *args, **options):
        """Entry point for command"""
        # Print message to the screen as the command is first executed.
        self.stdout.write("Waiting for database...")
        # Boolean value to check if the database is ready yet.
        db_up = False
        while db_up is False:
            try:
                """If the database is not ready, it throws an error
                and the except block will run."""
                self.check(databases=["default"])
                # If database is ready, the while loop is done.
                db_up = True
            except (Psycopg2OpError, OperationalError):
                # Print to screen when database is not ready yet.
                self.stdout.write("Database unavailable, waiting 1 second...")
                # Wait then run while loop again to check if database is ready.
                time.sleep(1)

        # After while loop is done, print success message.
        self.stdout.write(self.style.SUCCESS("Database available!"))
