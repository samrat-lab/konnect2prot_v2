# Konnect2Prot 2.0 - Web Application

## Setup and Run Instructions

### Prerequisites
- Ensure you have Python (>=3.8) installed.
- Install virtual environment support if not already available:
  ```bash
  pip install virtualenv
  ```

### Installation Steps

1. **Create and activate a virtual environment:**
   ```bash
   python -m venv venv
   source venv/bin/activate  # On macOS/Linux
   venv\Scripts\activate  # On Windows
   ```

2. **Install required dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

3. **Run the Flask application:**
   ```bash
   python app.py
   ```

### File Structure Overview
- `app.py`: Main entry point of the application.
- `requirements.txt`: Lists required Python packages.
- `templates/`: Contains HTML templates for rendering web pages.
- `static/`: Stores static files (CSS, JavaScript, images).
- `pages/`: Handles different sections of the app.
- `utils/`: Contains helper functions.
- `assets/`: Stores additional resources.

### Accessing the Application
Once the app is running, open a web browser and navigate to:
```
http://ip_address
```

### Deactivating the Virtual Environment
After usage, deactivate the virtual environment:
```bash
deactivate
```

This completes the setup process for running Konnect2Prot 2.0 locally.

