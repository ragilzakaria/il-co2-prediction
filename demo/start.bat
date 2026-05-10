@echo off
REM CO2 Solubility Prediction Demo - Startup Script (Windows)

echo 🚀 Starting CO₂ Solubility Prediction Demo...
echo.

REM Check if virtual environment exists
if not exist "venv" (
    echo 📦 Creating virtual environment...
    python -m venv venv
)

REM Activate virtual environment
echo 🔧 Activating virtual environment...
call venv\Scripts\activate.bat

REM Install dependencies if needed
if not exist "venv\.installed" (
    echo 📥 Installing dependencies...
    pip install -r requirements.txt
    type nul > venv\.installed
)

REM Check if model files exist
if not exist "..\models\best_random_forest_model.pkl" (
    echo ❌ Error: Model file not found!
    echo    Expected: ..\models\best_random_forest_model.pkl
    pause
    exit /b 1
)

REM Start the server
echo.
echo ✅ Starting Flask server...
echo 🌐 Open http://localhost:5000 in your browser
echo.
python app.py

pause

