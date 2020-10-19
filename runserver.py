"""
This script runs the FlaskWebProject1 application using a development server.
"""

from os import environ
import sys
sys.path.append('./')

from web_app import app
from web_app import views
from web_app.api import to
# from flask_script import Manager

# manager=Manager(app)

if __name__ == '__main__':

    HOST = environ.get('SERVER_HOST', '192.168.110.227')
    try:
        PORT = int(environ.get('SERVER_PORT', '5555'))
    except ValueError:
        PORT = 5555
    app.run(HOST, PORT)
