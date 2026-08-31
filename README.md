
# Chemistry Calculator (v1.0.0)
A zero-dependency Python tool and batch processor for chemistry and materials science. It solves ideal gas laws, molarity, dilution, liquid density, and stoichiometry equations with automatic unit conversions and built-in chemical formula parsing.

## Key Features
  It runs entirely on the standard Python 3 library so you don't have to worry about pip installs or virtual environments. 
  
  You can run it as an interactive menu for quick calculations or pass a CSV file to process hundreds of equations instantly. 
  
  The script parses chemical formulas like H2O and CuSO4 on the fly so you rarely have to look up molar masses manually. 
 
  It normalizes your units automatically allowing you to mix and match inputs like grams, atmospheres, and Celsius without doing the conversion math yourself. 
  
  The tool keeps a running history of your session so you can easily review your previous calculations. 

## Usage
  Just download the script and run it directly in your terminal.
  
  To start the interactive menu, run python3 gas_mole_calculator.py and follow the prompts.  
  
  If you want to use the batch processor, you can run python3 gas_mole_calculator.py --batch your_file.csv to process a spreadsheet of calculations. You can easily generate a sample     CSV template by selecting option 6 from the main menu.  

## Limitations in Version 1.0
  I made a few intentional trade-offs in this first version to keep everything contained in a single lightweight file. 
  
  The formula parser handles basic structures fine but it will throw an error if you try to use parentheses for groups or hydrate notation. 
  
  The gas solver strictly uses the ideal gas law and won't give accurate results for real gases at extreme pressures or very low temperatures. 
  
  The atomic masses are hardcoded to standard averages and won't accommodate specific isotopic calculations. 
  
  The code is structured as a standalone CLI tool rather than a modular package you can import into other Python apps. 

## Testing
  The easiest way to verify that the calculator and unit conversions are working correctly on your machine is to use the batch engine. Generate the sample CSV file from the main menu    and then run it through the batch processor. If it finishes without throwing any exceptions and the outputs look right, the core logic is fully functional.

## Contributing
  If you'd like to help add support for nested formulas, real gas equations, or a proper test suite, feel free to open a pull request. The only strict rule for this project is that we   don't add any external dependencies.

Why I Built This
  I created this tool to fix a constant annoyance during lab work and chemistry assignments: having to switch between a calculator, periodic table tabs, and unit conversion sites just   to work through routine math. Most software solutions were either clunky web tools or required installing heavy scientific libraries for basic tasks. I wanted a single script that     is fast, completely self-contained, and flexible enough to handle quick interactive checks or process a whole batch of lab data in seconds.

<img width="456" height="368" alt="Screenshot 2026-08-31 at 8 59 28 AM" src="https://github.com/user-attachments/assets/fbecaf71-61f3-4a5b-8707-8f226881c745" />
