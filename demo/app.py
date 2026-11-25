"""
Flask server for CO2 Solubility Prediction Demo
Serves the trained Random Forest model via REST API
"""

from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS
import joblib
import numpy as np
import pandas as pd
import os
import warnings

app = Flask(__name__, 
            static_folder='static',
            static_url_path='/static')
CORS(app)  # Enable CORS for local development

# Load model, scaler, and feature names
MODEL_PATH = '../models/thesis/best_random_forest_model.pkl'
SCALER_PATH = '../models/thesis/scaler.pkl'
FEATURE_NAMES_PATH = '../models/thesis/feature_names.pkl'

model = None
scaler = None
feature_names = None

def load_model():
    """Load the trained model, scaler, and feature names"""
    global model, scaler, feature_names
    
    try:
        # Suppress version warnings for pickle compatibility
        warnings.filterwarnings('ignore', category=UserWarning, module='sklearn')
        
        # Check if files exist
        if not os.path.exists(MODEL_PATH):
            print(f"❌ Model file not found: {MODEL_PATH}")
            print(f"   Current directory: {os.getcwd()}")
            return False
        if not os.path.exists(SCALER_PATH):
            print(f"❌ Scaler file not found: {SCALER_PATH}")
            return False
        if not os.path.exists(FEATURE_NAMES_PATH):
            print(f"❌ Feature names file not found: {FEATURE_NAMES_PATH}")
            return False
        
        model = joblib.load(MODEL_PATH)
        scaler = joblib.load(SCALER_PATH)
        feature_names = joblib.load(FEATURE_NAMES_PATH)
        print(f"✅ Model loaded successfully. Expects {len(feature_names)} features.")
        return True
    except Exception as e:
        print(f"❌ Error loading model: {e}")
        import traceback
        traceback.print_exc()
        return False

@app.route('/')
def index():
    """Serve the main HTML page"""
    return send_from_directory('static', 'index.html')

@app.route('/port')
def get_port():
    """Return the port number the server is running on"""
    import socket
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.connect(('127.0.0.1', 0))
    port = sock.getsockname()[1]
    sock.close()
    # Actually, let's just return the port from the request
    return jsonify({'port': request.environ.get('SERVER_PORT', 'unknown')})

@app.route('/predict', methods=['POST'])
def predict():
    """Make a prediction using the loaded model"""
    if model is None or scaler is None or feature_names is None:
        return jsonify({'error': 'Model not loaded'}), 500
    
    try:
        data = request.json
        
        # Extract input parameters
        temp_k = float(data.get('Temperature_K', 298.15))
        pressure_kpa = float(data.get('Pressure_kPa', 1000))
        
        # Cation properties (defaults for [BMIM])
        cation_energy = float(data.get('cation_Energy_Hartree', -423.0426))
        cation_homo = float(data.get('cation_HOMO_Eh', -0.416496))
        cation_lumo = float(data.get('cation_LUMO_Eh', -0.216006))
        cation_gap = float(data.get('cation_HOMO_LUMO_Gap_Eh', 0.20049))
        cation_dipole = float(data.get('cation_Dipole_Moment_Debye', 1.9587))
        
        # Anion properties
        anion_energy = float(data.get('anion_Energy_Hartree', -856.0656))
        anion_homo = float(data.get('anion_HOMO_Eh', 0.003031))
        anion_lumo = float(data.get('anion_LUMO_Eh', 0.053606))
        anion_gap = float(data.get('anion_HOMO_LUMO_Gap_Eh', 0.050575))
        anion_dipole = float(data.get('anion_Dipole_Moment_Debye', 13.4882))
        
        # Atom counts (defaults for [BMIM] cation)
        cation_atoms = data.get('cation_atoms', {'C': 8, 'H': 15, 'N': 2})
        anion_atoms = data.get('anion_atoms', {'C': 18, 'H': 33, 'O': 2})
        
        # Build feature vector
        feature_dict = {}
        
        # Initialize all features to 0
        for feat in feature_names:
            feature_dict[feat] = 0.0
        
        # Set conditions
        feature_dict['Temperature_K'] = temp_k
        feature_dict['Pressure_kPa'] = pressure_kpa
        
        # Set cation descriptors
        feature_dict['cation_Energy_Hartree'] = cation_energy
        feature_dict['cation_HOMO_Eh'] = cation_homo
        feature_dict['cation_LUMO_Eh'] = cation_lumo
        feature_dict['cation_HOMO_LUMO_Gap_Eh'] = cation_gap
        feature_dict['cation_Dipole_Moment_Debye'] = cation_dipole
        
        # Set anion descriptors
        feature_dict['anion_Energy_Hartree'] = anion_energy
        feature_dict['anion_HOMO_Eh'] = anion_homo
        feature_dict['anion_LUMO_Eh'] = anion_lumo
        feature_dict['anion_HOMO_LUMO_Gap_Eh'] = anion_gap
        feature_dict['anion_Dipole_Moment_Debye'] = anion_dipole
        
        # Set atom counts
        for element, count in cation_atoms.items():
            key = f'cation_num_{element}'
            if key in feature_dict:
                feature_dict[key] = float(count)
        
        for element, count in anion_atoms.items():
            key = f'anion_num_{element}'
            if key in feature_dict:
                feature_dict[key] = float(count)
        
        # Create DataFrame with exact feature order
        X_pred = pd.DataFrame([feature_dict])[feature_names]
        
        # Scale and predict
        X_pred_scaled = scaler.transform(X_pred)
        prediction = model.predict(X_pred_scaled)[0]
        
        return jsonify({
            'prediction': float(prediction),
            'temperature_K': temp_k,
            'pressure_kPa': pressure_kpa,
            'status': 'success'
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 400

@app.route('/health', methods=['GET'])
def health():
    """Health check endpoint"""
    if model is None:
        return jsonify({'status': 'not_ready', 'message': 'Model not loaded'}), 503
    return jsonify({
        'status': 'ready',
        'model_type': 'RandomForestRegressor',
        'n_features': len(feature_names) if feature_names else 0
    })

if __name__ == '__main__':
    import socket
    
    print("🚀 Starting CO2 Solubility Prediction Demo Server...")
    print(f"📁 Working directory: {os.getcwd()}")
    
    if load_model():
        # Find an available port starting from 5000
        port = 5000
        max_attempts = 10
        
        for attempt in range(max_attempts):
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            sock.settimeout(1)
            result = sock.connect_ex(('127.0.0.1', port))
            sock.close()
            
            if result != 0:  # Port is available
                break
            else:
                if attempt == 0:
                    print(f"⚠️  Port {port} is in use. Trying alternative ports...")
                port += 1
        
        if port >= 5000 + max_attempts:
            print(f"❌ Could not find an available port. Please free up ports 5000-{5000+max_attempts-1}")
            print("💡 On macOS, try disabling 'AirPlay Receiver' in System Settings > General > AirDrop & Handoff")
            exit(1)
        
        print(f"\n{'='*60}")
        print(f"✅ Server is ready!")
        print(f"🌐 URL: http://localhost:{port}")
        print(f"🌐 Alternative: http://127.0.0.1:{port} (use if Private Relay blocks localhost)")
        print(f"📖 Open one of these URLs in your browser")
        print(f"{'='*60}\n")
        app.run(debug=True, host='0.0.0.0', port=port, use_reloader=False)
    else:
        print("❌ Failed to load model. Please check file paths.")

