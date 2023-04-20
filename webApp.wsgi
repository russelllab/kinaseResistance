import os, sys
sys.path.insert(0, "/var/www/flask_apps/kinaseResistance/venv/lib/python3.6/site-packages/")

print (os.getcwd())
from webApp.routes import app as application
os.chdir('webApp')
