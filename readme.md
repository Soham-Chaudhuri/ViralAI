# ViralAI - Virus & Bacteria Gene Sequence Analyzer and Predictor

ViralAI is an advanced web application designed to analyze viral and bacterial gene sequences, predict the type of organism, and visualize protein structures. This tool leverages machine learning and bioinformatics to provide researchers with insights into genomes, enabling detailed sequence analysis and predictive modeling.

## 🚀 Features
- **Organism Prediction**: Predicts whether a sequence belongs to a virus or bacteria based on trained models.
- **Protein Structure Prediction**: Displays the predicted 3D protein structure using Py3DMol.
- **Sequence Analysis**: 
  - **GC vs AT Content Plot**: Visualize nucleotide composition.
  - **Hydrophobic vs Hydrophilic Score Plot**: Analyze amino acid properties.
  - **Molecular Weight vs Sequence Entropy Plot**: Explore relationships between molecular properties and sequence complexity.
- **Interactive Visualizations**: Generate real-time plots and 3D models in the browser.

## 🎥 Demo
Check out our video demonstration on YouTube: [ViralAI Demo](https://youtu.be/QRQzjd0yQOg)

## 🛠️ Technologies Used
### Frontend:
- **ReactJS**: Dynamic and interactive user interface.
- **Py3DMol**: 3D protein structure visualization in the browser.
- **Chart.js**: Plotting GC/AT content, hydrophobicity, and other metrics.

### Backend:
- **FastAPI**: High-performance Python backend for API development.
- **Python**: Core language for sequence analysis and ML integration.
- **Scikit-learn**: Trained ML models for classification.
- **Matplotlib**: For generating high-quality visualizations.
- **Biopython**: Sequence parsing and biological analysis.

## 🧪 Usage

1. **Input a Gene Sequence**:
   Paste your viral or bacterial gene sequence into the input field.

2. **Prediction**:
   The machine learning model will classify the sequence as one of the following:

   ```json
   {
       0: "Bacillus Subtilis",
       1: "Chikungunya Virus",
       2: "Clostridium Botulinum",
       3: "Escherichia Coli",
       4: "Herpes Simplex Virus 1",
       5: "Human Mastadenovirus C",
       6: "Human Papillomavirus",
       7: "Influenza A Virus",
       8: "Listeria Monocytogenes",
       9: "Pseudomonas Aeruginosa",
       10: "Staphylococcus Aureus",
       11: "Zika Virus"
   }

3. **Sequence Analysis**:
    - GC vs AT Content Plot: Visualize the ratio of GC content vs AT content in the input sequence.
    - Hydrophobic vs Hydrophilic Score Plot: Analyze the hydrophobicity and hydrophilicity of the amino acid sequence.
    - Molecular Weight vs Sequence Entropy Plot: Explore the relationship between molecular weight and the sequence's entropy.

4. **Protein Structure Visualization**:
    - A predicted 3D protein structure will be rendered in the browser using Py3DMol.
