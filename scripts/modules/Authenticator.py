import google.auth.transport.requests
import google.auth

from time import sleep
import os


# Needs pip library `requests`
class GoogleToken:
    def __init__(self):
        self.creds, self.project = google.auth.default()
        self.auth_req = google.auth.transport.requests.Request()

    def get_token(self):
        if (not self.creds.valid) or self.creds.expired:
            print("Refreshing token...")
            self.creds.refresh(self.auth_req)

        return self.creds.token

    def update_environment(self):
        if (not self.creds.valid) or self.creds.expired:
            print("Updating environment token...")
            self.creds.refresh(self.auth_req)

            os.environ["GCS_OAUTH_TOKEN"] = self.creds.token

    def test_expiration(self):
        for i in range(15):
            print(i)
            self.update_environment()
            sleep(60*10)
