# ASMRDE: Adaptive Social Mobility-Restructuring Differential Evolution

A MATLAB implementation of the Adaptive Social Mobility-Restructuring Differential Evolution algorithm for global optimization problems.

## Paper Information

**Title:** Adaptive social mobility-restructuring differential evolution for global optimization  
**Paper Link:** https://authors.elsevier.com/a/1li3x3PiGTXL4A

## Quick Start

### Prerequisites
- MATLAB (recommended version R2018a or later)
- Parallel Computing Toolbox (optional, for parallel execution)

### Running the Algorithm

1. Clone or download this repository
2. Open MATLAB and navigate to the ASMRDE folder
3. Run the main script:
   ```matlab
   CEC2017RUN
   ```

### Basic Configuration

By default, the algorithm runs on CEC2017 benchmark with 10 dimensions. To change the dimension, modify the `SearchDimension` parameter in `CEC2017RUN.m`:

```matlab
SearchDimension = 10;  % Change to desired dimension (10, 30, 50, or 100)
```

## File Structure

```
ASMRDE/
念岸岸 ASMRDE.m           # Main algorithm implementation
念岸岸 CEC2017RUN.m       # Algorithm runner and benchmark testing
念岸岸 cec17_func.mexw64  # CEC2017 benchmark functions (compiled)
念岸岸 input_data/        # CEC2017 test data files
岫   念岸岸 M_1_D10.txt
岫   念岸岸 M_1_D100.txt
岫   弩岸岸 ...
弩岸岸 README.md          # This file
```

### Key Files Description

- **`ASMRDE.m`**: Core algorithm implementation containing the Adaptive Social Mobility-Restructuring Differential Evolution logic
- **`CEC2017RUN.m`**: Main execution script that handles parallel computing, parameter settings, and result collection
- **`cec17_func.mexw64`**: Compiled CEC2017 benchmark functions
- **`input_data/`**: Contains transformation matrices and shift data for CEC2017 test functions

## Customization

### Using Different Benchmark Sets

To use other CEC benchmark sets (e.g., CEC2013, CEC2014):

1. Replace the benchmark files:
   - Update `cec17_func.mexw64` with the corresponding benchmark file (e.g., `cec13_func.mexw64`)
   - Replace the `input_data/` folder with the appropriate test data

2. Modify the function call in `ASMRDE.m`:
   ```matlab
   % Line 36 and 165: Change cec17_func to your benchmark function
   POPS(1:PopSize,SearchDimension+1) = cec13_func(POPS(1:PopSize,1:SearchDimension)',FuncNo)';
   U(1:PopSize,SearchDimension+1) = cec13_func(U(1:PopSize,1:SearchDimension)',FuncNo);
   ```

### Algorithm Parameters

Key parameters can be adjusted in `ASMRDE.m`:

- **Population sizes** (lines 21-28 in `CEC2017RUN.m`):
  - 10D: 90 individuals
  - 30D: 150 individuals  
  - 50D: 210 individuals
  - 100D: 260 individuals

- **Loop period** (line 72): Social restructuring frequency
- **Diversity threshold** (line 75): Population diversity control
- **Temperature parameters** (line 86): Boltzmann distribution parameter

## Output

The algorithm generates:
- **`.mat` files**: Complete results with detailed statistics
- **`.xlsx` files**: Mean values for each test function
- **Console output**: Real-time progress and statistics

Results include:
- Mean, Best, Worst fitness values
- Standard deviation
- Convergence curves
- Population diversity metrics

## Contributing

If you encounter any issues or have suggestions for improvement:

1. **Preferred**: Open an issue on GitHub
2. **Alternative**: Contact via email

## Citation

If you use this code in your research, please cite the original paper:

```bibtex
@article{ASMRDE2024,
  title={Adaptive social mobility-restructuring differential evolution for global optimization},
  author={[Yiwen Zhuo]},
  journal={[Expert Systems with Applications]},
  year={2026},
  url={https://authors.elsevier.com/a/1li3x3PiGTXL4A}
}
```

## System Requirements

- **MATLAB**: R2018a or later recommended
- **Memory**: At least 4GB RAM for larger dimensions
- **Storage**: ~50MB for CEC2017 data files
- **Optional**: Parallel Computing Toolbox for faster execution

---


For more details about the algorithm mechanics and experimental results, please refer to the original paper. 
