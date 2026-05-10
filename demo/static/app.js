// API endpoint - will be detected from the page URL
const API_URL = window.location.origin;

// DOM elements
const form = document.getElementById('predictionForm');
const predictBtn = document.getElementById('predictBtn');
const resetBtn = document.getElementById('resetBtn');
const resultCard = document.getElementById('resultCard');

// Default values
const defaults = {
    Temperature_K: 298.15,
    Pressure_kPa: 1000,
    cation_Energy_Hartree: -423.0426,
    cation_HOMO_Eh: -0.416496,
    cation_LUMO_Eh: -0.216006,
    cation_HOMO_LUMO_Gap_Eh: 0.20049,
    cation_Dipole_Moment_Debye: 1.9587,
    anion_Energy_Hartree: -856.0656,
    anion_HOMO_Eh: 0.003031,
    anion_LUMO_Eh: 0.053606,
    anion_HOMO_LUMO_Gap_Eh: 0.050575,
    anion_Dipole_Moment_Debye: 13.4882
};

// Default atom counts (empty - user adds manually)
const defaultCationAtoms = {};
const defaultAnionAtoms = {};

// Track selected atoms
let cationSelectedAtoms = new Set();
let anionSelectedAtoms = new Set();

// Initialize form with defaults
function initializeForm() {
    Object.keys(defaults).forEach(key => {
        const input = document.getElementById(key);
        if (input) {
            input.value = defaults[key];
        }
    });
    
    // Initialize atom inputs with defaults
    initializeAtomInputs('cation', defaultCationAtoms);
    initializeAtomInputs('anion', defaultAnionAtoms);
}

// Initialize atom inputs
function initializeAtomInputs(type, atoms) {
    const container = document.getElementById(`${type}_atom_inputs`);
    container.innerHTML = '';
    
    if (type === 'cation') {
        cationSelectedAtoms = new Set(Object.keys(atoms));
    } else {
        anionSelectedAtoms = new Set(Object.keys(atoms));
    }
    
    Object.keys(atoms).forEach(atom => {
        addAtomInput(type, atom, atoms[atom]);
    });
    
    // Update dropdown after initialization
    updateAtomDropdown(type);
}

// Add atom input field
function addAtomInput(type, atom, value = 0) {
    const container = document.getElementById(`${type}_atom_inputs`);
    const selectedSet = type === 'cation' ? cationSelectedAtoms : anionSelectedAtoms;
    const select = document.getElementById(`${type}_atom_select`);
    
    if (!atom || atom === '') {
        return; // No atom selected
    }
    
    if (selectedSet.has(atom)) {
        // Atom already added - show feedback
        const option = Array.from(select.options).find(opt => opt.value === atom);
        if (option) {
            select.value = atom; // Highlight it
            setTimeout(() => {
                select.value = '';
            }, 500);
        }
        return; // Already added
    }
    
    selectedSet.add(atom);
    
    const item = document.createElement('div');
    item.className = 'atom-input-item';
    item.id = `${type}_atom_${atom}_container`;
    item.innerHTML = `
        <label>${atom}:</label>
        <input type="number" id="${type}_${atom}" value="${value}" min="0" step="1" class="atom-count-input">
        <button type="button" class="btn-remove-atom" onclick="removeAtom('${type}', '${atom}')" title="Remove ${atom}">×</button>
    `;
    
    container.appendChild(item);
    
    // Remove from dropdown or disable it
    updateAtomDropdown(type);
    
    // Reset dropdown selection
    select.value = '';
}

// Update dropdown to disable/hide already selected atoms
function updateAtomDropdown(type) {
    const select = document.getElementById(`${type}_atom_select`);
    const selectedSet = type === 'cation' ? cationSelectedAtoms : anionSelectedAtoms;
    
    Array.from(select.options).forEach(option => {
        if (option.value && selectedSet.has(option.value)) {
            option.disabled = true;
            option.textContent = option.textContent + ' (already added)';
        } else if (option.value) {
            // Reset if atom was removed
            const originalText = option.textContent.replace(' (already added)', '');
            option.disabled = false;
            option.textContent = originalText;
        }
    });
}

// Remove atom input
function removeAtom(type, atom) {
    const container = document.getElementById(`${type}_atom_${atom}_container`);
    if (container) {
        container.remove();
        const selectedSet = type === 'cation' ? cationSelectedAtoms : anionSelectedAtoms;
        selectedSet.delete(atom);
        // Re-enable in dropdown
        updateAtomDropdown(type);
    }
}

// Setup atom selection handlers
function setupAtomHandlers() {
    // Cation atom handler
    document.getElementById('add_cation_atom').addEventListener('click', () => {
        const select = document.getElementById('cation_atom_select');
        const atom = select.value;
        if (atom && !select.options[select.selectedIndex].disabled) {
            addAtomInput('cation', atom, 0);
        } else if (atom) {
            alert(`${atom} is already added. Remove it first if you want to change the count.`);
        }
    });
    
    // Anion atom handler
    document.getElementById('add_anion_atom').addEventListener('click', () => {
        const select = document.getElementById('anion_atom_select');
        const atom = select.value;
        if (atom && !select.options[select.selectedIndex].disabled) {
            addAtomInput('anion', atom, 0);
        } else if (atom) {
            alert(`${atom} is already added. Remove it first if you want to change the count.`);
        }
    });
    
    // Allow Enter key to add atom
    document.getElementById('cation_atom_select').addEventListener('keypress', (e) => {
        if (e.key === 'Enter') {
            e.preventDefault();
            document.getElementById('add_cation_atom').click();
        }
    });
    
    document.getElementById('anion_atom_select').addEventListener('keypress', (e) => {
        if (e.key === 'Enter') {
            e.preventDefault();
            document.getElementById('add_anion_atom').click();
        }
    });
    
    // Update dropdowns on initialization
    updateAtomDropdown('cation');
    updateAtomDropdown('anion');
}

// Reset form to defaults
resetBtn.addEventListener('click', () => {
    initializeForm();
    resultCard.innerHTML = `
        <div class="result-placeholder">
            <p>👆 Enter parameters and click "Predict" to see results</p>
        </div>
    `;
});

// Handle form submission
form.addEventListener('submit', async (e) => {
    e.preventDefault();
    
    // Disable button and show loading
    predictBtn.disabled = true;
    predictBtn.innerHTML = '<span class="loading"></span> Predicting...';
    
    // Collect form data
    const formData = {
        Temperature_K: parseFloat(document.getElementById('Temperature_K').value),
        Pressure_kPa: parseFloat(document.getElementById('Pressure_kPa').value),
        cation_Energy_Hartree: parseFloat(document.getElementById('cation_Energy_Hartree').value),
        cation_HOMO_Eh: parseFloat(document.getElementById('cation_HOMO_Eh').value),
        cation_LUMO_Eh: parseFloat(document.getElementById('cation_LUMO_Eh').value),
        cation_HOMO_LUMO_Gap_Eh: parseFloat(document.getElementById('cation_HOMO_LUMO_Gap_Eh').value),
        cation_Dipole_Moment_Debye: parseFloat(document.getElementById('cation_Dipole_Moment_Debye').value),
        anion_Energy_Hartree: parseFloat(document.getElementById('anion_Energy_Hartree').value),
        anion_HOMO_Eh: parseFloat(document.getElementById('anion_HOMO_Eh').value),
        anion_LUMO_Eh: parseFloat(document.getElementById('anion_LUMO_Eh').value),
        anion_HOMO_LUMO_Gap_Eh: parseFloat(document.getElementById('anion_HOMO_LUMO_Gap_Eh').value),
        anion_Dipole_Moment_Debye: parseFloat(document.getElementById('anion_Dipole_Moment_Debye').value),
        cation_atoms: collectAtomCounts('cation'),
        anion_atoms: collectAtomCounts('anion')
    };
    
    try {
        // Make API request
        console.log('Making prediction request to:', `${API_URL}/predict`);
        console.log('Request data:', formData);
        
        const response = await fetch(`${API_URL}/predict`, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify(formData)
        });
        
        console.log('Response status:', response.status);
        const data = await response.json();
        console.log('Response data:', data);
        
        if (response.ok && data.status === 'success') {
            // Display result
            console.log('Displaying result:', data);
            displayResult(data);
        } else {
            // Display error
            console.error('Prediction failed:', data);
            displayError(data.error || 'Prediction failed');
        }
    } catch (error) {
        console.error('Error:', error);
        displayError(`Connection error: ${error.message}. Make sure the server is running on ${API_URL}`);
    } finally {
        // Re-enable button
        predictBtn.disabled = false;
        predictBtn.innerHTML = '🧪 Predict CO₂ Solubility';
    }
});

// Display prediction result
function displayResult(data) {
    console.log('displayResult called with:', data);
    
    if (!resultCard) {
        console.error('resultCard element not found!');
        return;
    }
    
    const prediction = data.prediction;
    const temp = data.temperature_K;
    const pressure = data.pressure_kPa;
    
    if (prediction === undefined || prediction === null) {
        console.error('Invalid prediction value:', prediction);
        displayError('Invalid prediction result received from server');
        return;
    }
    
    resultCard.innerHTML = `
        <div class="result-content">
            <div class="prediction-label">Predicted CO₂ Solubility</div>
            <div class="prediction-value">${prediction.toFixed(4)}</div>
            <div class="prediction-details">
                <div>Temperature: ${temp} K (${(temp - 273.15).toFixed(1)} °C)</div>
                <div>Pressure: ${pressure} kPa (${(pressure / 100).toFixed(1)} bar)</div>
                <div style="margin-top: 10px; font-size: 0.85rem; opacity: 0.8;">
                    Mole fraction of CO₂ in ionic liquid
                </div>
            </div>
        </div>
    `;
    
    // Scroll to result
    resultCard.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
    
    console.log('Result displayed successfully');
}

// Display error message
function displayError(message) {
    resultCard.innerHTML = `
        <div class="error-message">
            <strong>⚠️ Error</strong><br>
            ${message}
        </div>
    `;
}

// Check server health on load
async function checkServerHealth() {
    try {
        const response = await fetch(`${API_URL}/health`);
        const data = await response.json();
        
        if (data.status === 'ready') {
            console.log('✅ Server is ready');
        } else {
            console.warn('⚠️ Server not ready:', data.message);
        }
    } catch (error) {
        console.error('❌ Cannot connect to server:', error);
        resultCard.innerHTML = `
            <div class="error-message">
                <strong>⚠️ Server Not Available</strong><br>
                Please make sure the Flask server is running.<br>
                Run: <code>python app.py</code> in the demo directory
            </div>
        `;
    }
}

// Collect atom counts from inputs
function collectAtomCounts(type) {
    const atoms = {};
    const selectedSet = type === 'cation' ? cationSelectedAtoms : anionSelectedAtoms;
    
    selectedSet.forEach(atom => {
        const input = document.getElementById(`${type}_${atom}`);
        if (input) {
            const value = parseInt(input.value) || 0;
            if (value > 0) {
                atoms[atom] = value;
            }
        }
    });
    
    return atoms;
}

// Make removeAtom available globally
window.removeAtom = removeAtom;

// Initialize on page load
document.addEventListener('DOMContentLoaded', () => {
    console.log('Page loaded, initializing...');
    console.log('API_URL:', API_URL);
    console.log('resultCard element:', resultCard);
    
    if (!resultCard) {
        console.error('❌ resultCard element not found in DOM!');
        alert('Error: Result card element not found. Please refresh the page.');
        return;
    }
    
    setupAtomHandlers();
    initializeForm();
    checkServerHealth();
    
    console.log('Initialization complete');
});

