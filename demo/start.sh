#!/bin/bash

# CO2 Solubility Prediction Demo - Startup Script

echo "🚀 Starting CO₂ Solubility Prediction Demo..."
echo ""

# Check if virtual environment exists
if [ ! -d "venv" ]; then
    echo "📦 Creating virtual environment..."
    python3 -m venv venv
fi

# Activate virtual environment
echo "🔧 Activating virtual environment..."
source venv/bin/activate

# Install dependencies if needed
if [ ! -f "venv/.installed" ]; then
    echo "📥 Installing dependencies..."
    pip install --upgrade pip
    pip install -r requirements.txt
    if [ $? -eq 0 ]; then
        touch venv/.installed
    else
        echo "⚠️  Installation had issues. Trying with system packages..."
        echo "💡 Tip: If scikit-learn fails to install, try using Python 3.11 or 3.12"
        exit 1
    fi
fi

# Check if model files exist
if [ ! -f "../models/best_random_forest_model.pkl" ]; then
    echo "❌ Error: Model file not found!"
    echo "   Expected: ../models/best_random_forest_model.pkl"
    exit 1
fi

# Start the server
echo ""
echo "✅ Starting Flask server..."
echo "🌐 Open http://localhost:5000 in your browser"
echo ""
python app.py

