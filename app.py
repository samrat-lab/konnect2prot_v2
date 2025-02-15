import dash
import base64
import io
import pandas as pd
import dash_uploader as du
from dash import dcc, html
from flask import session, Flask, render_template ,request
from pymongo import MongoClient
#from flask_session import Session
from call_backs import register_callbacks
from dash.dependencies import Input, Output, State
from index import index_layout #layout as index_layout
import uuid 
import arrow
server = Flask(__name__)
client = MongoClient("mongodb://127.0.0.1:27017")
db = client.ppi
@server.route('/')
def index():
    return render_template('home_page.html')


@server.route('/feedback', methods=["GET", "POST"])
def form_entry():
    data = {}
    if request.method == "POST":
        data['Name'] = request.form['name']
        
        
        data['Email'] = request.form['email']
        data['Subject'] = request.form['subject']
        data['Message'] = request.form['message']
        
        data['Date_Time'] = arrow.now("Asia/Calcutta").format()
       
        db.feedback.insert_one(data)

    return render_template("home_page.html")

@server.route('/tutorial')
def tutorial():
    return render_template('k2p_tutorial.html')

app = dash.Dash(__name__, server=server, external_scripts=['https://cdn.tailwindcss.com'],routes_pathname_prefix='/dashboard/', suppress_callback_exceptions=True)
app.title="k2p_v2"

# Configure upload folder
UPLOAD_FOLDER_ROOT = r"data_store"
du.configure_upload(app, UPLOAD_FOLDER_ROOT)

# Configure server-side session storage
app.server.config["SECRET_KEY"] = "your-secret-key"  # Replace with your secret key
app.server.config["SESSION_TYPE"] = "filesystem"  # Store session data in the filesystem
#Session(app.server)

# Generate a static UUID if it doesn't exist in the session
  # Generate and store in session

# Set the main layout
app.layout = index_layout(app)
register_callbacks(app)

if __name__ == '__main__':
    app.run_server(host='ip_address', debug=True)
