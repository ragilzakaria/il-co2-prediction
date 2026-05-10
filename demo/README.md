# CO₂ Solubility Prediction Demo

A local web-based demo application for predicting CO₂ solubility in ionic liquids using a trained Random Forest model.

## Features

- 🌐 **Web Interface**: Modern, responsive HTML/CSS/JavaScript frontend
- 🧪 **Real-time Predictions**: Interactive form to input molecular descriptors
- 📊 **Model Information**: Display of model details and prediction results
- 🎨 **Beautiful UI**: Gradient design with smooth animations

## Prerequisites

- Python 3.8 or higher
- pip (Python package manager)

## Installation

1. **Navigate to the demo directory:**
   ```bash
   cd demo
   ```

2. **Create a virtual environment (recommended):**
   ```bash
   python -m venv venv
   
   # On macOS/Linux:
   source venv/bin/activate
   
   # On Windows:
   venv\Scripts\activate
   ```

3. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

## Running the Demo

### Quick Start (Recommended)

**On macOS/Linux:**
```bash
./start.sh
```

**On Windows:**
```batch
start.bat
```

The startup script will:
- Create a virtual environment (if needed)
- Install dependencies automatically
- Start the Flask server

### Manual Start

1. **Navigate to the demo directory:**
   ```bash
   cd demo
   ```

2. **Activate virtual environment:**
   ```bash
   # macOS/Linux:
   source venv/bin/activate
   
   # Windows:
   venv\Scripts\activate
   ```

3. **Start the Flask server:**
   ```bash
   python app.py
   ```

   You should see:
   ```
   🚀 Starting CO2 Solubility Prediction Demo Server...
   ✅ Model loaded successfully. Expects 248 features.
   🌐 Server starting on http://localhost:5000
   📖 Open http://localhost:5000 in your browser
   ```

4. **Open your web browser:**
   Navigate to: `http://localhost:5000`

5. **Make predictions:**
   - Adjust the input parameters (temperature, pressure, molecular descriptors)
   - Click "🧪 Predict CO₂ Solubility"
   - View the predicted CO₂ solubility (mole fraction)

## Input Parameters

### Experimental Conditions
- **Temperature (K)**: Operating temperature (default: 298.15 K = 25°C)
- **Pressure (kPa)**: Operating pressure (default: 1000 kPa = 10 bar)

### Cation Properties
- Energy (Hartree)
- HOMO (Eh)
- LUMO (Eh)
- HOMO-LUMO Gap (Eh)
- Dipole Moment (Debye)

### Anion Properties
- Energy (Hartree)
- HOMO (Eh)
- LUMO (Eh)
- HOMO-LUMO Gap (Eh)
- Dipole Moment (Debye)

### Atom Counts
- Cation: C, H, N (and other elements as needed)
- Anion: C, H, O (and other elements as needed)

## Model Information

- **Model Type**: Random Forest Regressor
- **Number of Features**: 248 quantum descriptors
- **Target**: CO₂ solubility (mole fraction)
- **Training Data**: Experimental dataset with 11,256 samples

## API Endpoints

### `POST /predict`
Make a prediction with custom parameters.

**Request Body:**
```json
{
  "Temperature_K": 298.15,
  "Pressure_kPa": 1000,
  "cation_Energy_Hartree": -423.0426,
  "cation_HOMO_Eh": -0.416496,
  "cation_LUMO_Eh": -0.216006,
  "cation_HOMO_LUMO_Gap_Eh": 0.20049,
  "cation_Dipole_Moment_Debye": 1.9587,
  "anion_Energy_Hartree": -856.0656,
  "anion_HOMO_Eh": 0.003031,
  "anion_LUMO_Eh": 0.053606,
  "anion_HOMO_LUMO_Gap_Eh": 0.050575,
  "anion_Dipole_Moment_Debye": 13.4882,
  "cation_atoms": {"C": 8, "H": 15, "N": 2},
  "anion_atoms": {"C": 18, "H": 33, "O": 2}
}
```

**Response:**
```json
{
  "prediction": 0.1292,
  "temperature_K": 298.15,
  "pressure_kPa": 1000,
  "status": "success"
}
```

### `GET /health`
Check if the server and model are ready.

**Response:**
```json
{
  "status": "ready",
  "model_type": "RandomForestRegressor",
  "n_features": 248
}
```

## Troubleshooting

### scikit-learn installation fails (Python 3.13)
If you're using Python 3.13 and scikit-learn fails to install, the requirements.txt has been updated to use scikit-learn >=1.4.0 which has pre-built wheels for Python 3.13. If you still encounter issues:

1. **Option 1**: Use Python 3.11 or 3.12 (recommended for maximum compatibility)
   ```bash
   python3.11 -m venv venv
   # or
   python3.12 -m venv venv
   ```

2. **Option 2**: Install dependencies manually with latest versions
   ```bash
   pip install --upgrade pip setuptools wheel
   pip install flask flask-cors scikit-learn pandas numpy joblib
   ```

### Model files not found
Make sure the model files are in the correct location:
- `../models/best_random_forest_model.pkl`
- `../models/scaler.pkl`
- `../models/feature_names.pkl`

### Port already in use
The app automatically detects if port 5000 is in use and tries alternative ports (5001, 5002, etc.). If you see a port conflict:

1. **On macOS**: Disable 'AirPlay Receiver' in System Settings > General > AirDrop & Handoff
2. **Or**: The app will automatically use the next available port and display it in the terminal
3. **Or**: Manually specify a port by editing `app.py`:
   ```python
   app.run(debug=True, host='0.0.0.0', port=5001)  # Change to your preferred port
   ```

### CORS errors
The Flask server includes CORS support. If you encounter CORS issues, make sure `flask-cors` is installed.

### Connection Refused / Can't Connect to localhost

If you see "connection refused" or can't connect to the server:

1. **Check if Private Relay is enabled** (macOS):
   - Go to System Settings > [Your Name] > iCloud > Private Relay
   - **Temporarily disable Private Relay** or add an exception for localhost
   - Private Relay can interfere with localhost connections

2. **Try using `127.0.0.1` instead of `localhost`**:
   - Instead of: `http://localhost:5001`
   - Use: `http://127.0.0.1:5001`
   - This bypasses DNS resolution that Private Relay might interfere with

3. **Check the correct port**:
   - Look at the terminal output when starting the server
   - It will show: `🌐 Server starting on http://localhost:XXXX`
   - Use that exact port number

4. **Verify the server is running**:
   ```bash
   # Check if server is running
   curl http://127.0.0.1:5001/health
   # Should return: {"status":"ready",...}
   ```

5. **Check firewall settings**:
   - macOS: System Settings > Network > Firewall
   - Make sure Python/Flask is allowed through the firewall

## File Structure

```
demo/
├── app.py              # Flask server
├── requirements.txt    # Python dependencies
├── README.md           # This file
└── static/
    ├── index.html      # Main HTML page
    ├── styles.css      # CSS styling
    └── app.js          # JavaScript logic
```

## License

This demo is part of the CO₂ solubility prediction project.

