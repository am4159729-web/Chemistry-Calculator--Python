#!/usr/bin/env python3
"""
================================================================================
VERSATILE MASS, VOLUME, AND GAS MOLE CALCULATOR
================================================================================
"""
import sys
import os
import csv
import math
import re
from dataclasses import dataclass
from typing import Dict, Tuple, Optional, List, Any


# ==============================================================================
# PHYSICAL CONSTANTS & UNIT DEFINITIONS (SI Base Units Internal)
# ==============================================================================

# Universal Gas Constant R in SI units: J / (mol * K) == (Pa * m^3) / (mol * K)
R_SI = 8.31446261815324

# Built-in periodic table (standard atomic weights in g/mol).
# This supports the most commonly used elements in general chemistry and
# materials calculations without requiring an external lookup.
ATOMIC_MASSES: Dict[str, float] = {
    'H': 1.008, 'He': 4.0026, 'Li': 6.94, 'Be': 9.0122, 'B': 10.81,
    'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180,
    'Na': 22.989, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.085, 'P': 30.974,
    'S': 32.06, 'Cl': 35.45, 'Ar': 39.95, 'K': 39.098, 'Ca': 40.078,
    'Sc': 44.956, 'Ti': 47.867, 'V': 50.942, 'Cr': 51.996, 'Mn': 54.938,
    'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693, 'Cu': 63.546, 'Zn': 65.38,
    'Ga': 69.723, 'Ge': 72.630, 'As': 74.922, 'Se': 78.971, 'Br': 79.904,
    'Kr': 83.798, 'Rb': 85.468, 'Sr': 87.62, 'Y': 88.906, 'Zr': 91.224,
    'Nb': 92.906, 'Mo': 95.95, 'Ru': 101.07, 'Rh': 102.91, 'Pd': 106.42,
    'Ag': 107.87, 'Cd': 112.41, 'In': 114.82, 'Sn': 118.71, 'Sb': 121.76,
    'Te': 127.60, 'I': 126.90, 'Xe': 131.29, 'Cs': 132.91, 'Ba': 137.33,
    'W': 183.84, 'Pt': 195.08, 'Au': 196.97, 'Hg': 200.59, 'Pb': 207.2,
    'Bi': 208.98, 'U': 238.03
}

# Base Unit Mappings -> Factor to multiply by to convert unit to SI Base Unit
# Base units: Mass -> grams (g), Volume -> cubic meters (m^3), 
#             Moles -> moles (mol), Pressure -> Pascals (Pa)

MASS_UNITS: Dict[str, float] = {
    'g': 1.0, 'gram': 1.0, 'grams': 1.0,
    'kg': 1000.0, 'kilogram': 1000.0, 'kilograms': 1000.0,
    'mg': 0.001, 'milligram': 0.001, 'milligrams': 0.001,
    'ug': 1e-6, 'microgram': 1e-6, 'micrograms': 1e-6,
    'lb': 453.59237, 'lbs': 453.59237, 'pound': 453.59237, 'pounds': 453.59237,
    'oz': 28.349523125, 'ounce': 28.349523125, 'ounces': 28.349523125
}

VOLUME_UNITS: Dict[str, float] = {
    'm3': 1.0, 'm^3': 1.0, 'cubic_meter': 1.0, 'cubic_meters': 1.0,
    'l': 0.001, 'liter': 0.001, 'liters': 0.001, 'litre': 0.001, 'litres': 0.001,
    'ml': 1e-6, 'milliliter': 1e-6, 'milliliters': 1e-6, 'cc': 1e-6, 'cm3': 1e-6, 'cm^3': 1e-6,
    'gal': 0.003785411784, 'gallon': 0.003785411784, 'gallons': 0.003785411784
}

MOLE_UNITS: Dict[str, float] = {
    'mol': 1.0, 'mole': 1.0, 'moles': 1.0,
    'mmol': 0.001, 'millimole': 0.001, 'millimoles': 0.001,
    'kmol': 1000.0, 'kilomole': 1000.0, 'kilomoles': 1000.0,
    'umol': 1e-6, 'micromole': 1e-6, 'micromoles': 1e-6
}

PRESSURE_UNITS: Dict[str, float] = {
    'pa': 1.0, 'pascal': 1.0, 'pascals': 1.0,
    'kpa': 1000.0, 'kilopascal': 1000.0, 'kilopascals': 1000.0,
    'mpa': 1e6, 'megapascal': 1e6, 'megapascals': 1e6,
    'atm': 101325.0, 'atmosphere': 101325.0, 'atmospheres': 101325.0,
    'torr': 101325.0 / 760.0, 'torrs': 101325.0 / 760.0, 'mmhg': 101325.0 / 760.0,
    'bar': 100000.0, 'bars': 100000.0,
    'mbar': 100.0, 'millibar': 100.0,
    'psi': 6894.757293168
}

MOLARITY_UNITS: Dict[str, float] = {
    'm': 1.0, 'molar': 1.0, 'mol/l': 1.0, 'mol/dm3': 1.0,
    'mm': 0.001, 'millimolar': 0.001, 'mmol/l': 0.001
}

# Density factors convert to the internal base unit of g/mL.
DENSITY_UNITS: Dict[str, float] = {
    'g/ml': 1.0, 'g/cm3': 1.0, 'g/cc': 1.0, 'kg/l': 1.0,
    'kg/m3': 0.001, 'g/l': 0.001
}


@dataclass
class CalculationEntry:
    operation: str
    inputs: Dict[str, str]
    result: str
    explanation: str


# Global Session History
HISTORY: List[CalculationEntry] = []


# ==============================================================================
# CONVERSION ENGINE
# ==============================================================================

def normalize_unit_string(unit_str: str) -> str:
    """Normalize input unit strings (strip whitespace, lower-case)."""
    return str(unit_str).strip().lower()

def get_unit_factor(unit_map: Dict[str, float], unit: str, category: str) -> float:
    """Return a unit conversion factor with a clear, user-friendly error."""
    normalized = normalize_unit_string(unit)
    if normalized not in unit_map:
        supported = ", ".join(sorted(unit_map.keys()))
        raise ValueError(
            f"Unsupported {category} unit '{unit}'. "
            f"Supported units: {supported}"
        )
    return unit_map[normalized]

def calculate_molar_mass(input_value: str) -> float:
    """
    Parse a numeric molar mass or a simple chemical formula.

    Supported formula syntax includes element symbols with optional integer
    subscripts, for example H2O, CuSO4, and BaTiO3. Parentheses and hydrate
    notation are intentionally rejected rather than silently miscalculated.
    """
    formula = str(input_value).strip()
    if not formula:
        raise ValueError("Molar mass or chemical formula cannot be empty.")

    try:
        numeric_mass = float(formula)
    except ValueError:
        numeric_mass = None
    else:
        if not math.isfinite(numeric_mass) or numeric_mass <= 0:
            raise ValueError("Molar mass must be a finite value greater than zero.")
        return numeric_mass

    if not re.fullmatch(r"(?:[A-Z][a-z]*\d*)+", formula):
        raise ValueError(
            f"Could not parse '{formula}'. Enter a positive number or a "
            "formula such as H2O, CuSO4, or BaTiO3."
        )

    total_mass = 0.0
    for element, count_text in re.findall(r"([A-Z][a-z]*)(\d*)", formula):
        if element not in ATOMIC_MASSES:
            raise ValueError(
                f"Unknown element '{element}'. Check the chemical symbol."
            )
        count = int(count_text) if count_text else 1
        if count <= 0:
            raise ValueError("Element counts in a formula must be greater than zero.")
        total_mass += ATOMIC_MASSES[element] * count

    if total_mass <= 0:
        raise ValueError("Calculated molar mass must be greater than zero.")
    return total_mass

def parse_environment_input(value: str, unit: str, env_type: str) -> Tuple[float, str]:
    """Parse a numeric batch value or an STP/RT/SATP preset."""
    raw_value = str(value).strip().upper()
    if env_type == "temperature":
        if raw_value == "STP":
            return 0.0, "C"
        if raw_value in {"RT", "SATP"}:
            return 25.0, "C"
        try:
            return float(value), unit
        except ValueError:
            raise ValueError(
                f"Invalid temperature '{value}'. Use a number, STP, RT, or SATP."
            )
    if env_type == "pressure":
        if raw_value in {"STP", "RT", "SATP"}:
            return 1.0, "atm"
        try:
            return float(value), unit
        except ValueError:
            raise ValueError(
                f"Invalid pressure '{value}'. Use a number, STP, RT, or SATP."
            )
    raise ValueError("Environment type must be temperature or pressure.")

def convert_mass(value: float, from_unit: str, to_unit: str) -> float:
    grams = value * get_unit_factor(MASS_UNITS, from_unit, "mass")
    return grams / get_unit_factor(MASS_UNITS, to_unit, "mass")

def convert_volume(value: float, from_unit: str, to_unit: str) -> float:
    m3 = value * get_unit_factor(VOLUME_UNITS, from_unit, "volume")
    return m3 / get_unit_factor(VOLUME_UNITS, to_unit, "volume")

def convert_mole(value: float, from_unit: str, to_unit: str) -> float:
    moles = value * get_unit_factor(MOLE_UNITS, from_unit, "mole")
    return moles / get_unit_factor(MOLE_UNITS, to_unit, "mole")

def convert_pressure(value: float, from_unit: str, to_unit: str) -> float:
    pascals = value * get_unit_factor(PRESSURE_UNITS, from_unit, "pressure")
    return pascals / get_unit_factor(PRESSURE_UNITS, to_unit, "pressure")

def to_kelvin(value: float, unit: str) -> float:
    u = normalize_unit_string(unit)
    if u in ['k', 'kelvin']:
        return value
    elif u in ['c', 'celsius']:
        return value + 273.15
    elif u in ['f', 'fahrenheit']:
        return (value - 32.0) * (5.0 / 9.0) + 273.15
    else:
        raise ValueError(
            f"Unsupported temperature unit '{unit}'. "
            "Supported units: K, C, F."
        )

def from_kelvin(kelvin_val: float, target_unit: str) -> float:
    u = normalize_unit_string(target_unit)
    if u in ['k', 'kelvin']:
        return kelvin_val
    elif u in ['c', 'celsius']:
        return kelvin_val - 273.15
    elif u in ['f', 'fahrenheit']:
        return (kelvin_val - 273.15) * (9.0 / 5.0) + 32.0
    else:
        raise ValueError(
            f"Unsupported temperature unit '{target_unit}'. "
            "Supported units: K, C, F."
        )


# ==============================================================================
# CORE CORE CALCULATION LOGIC
# ==============================================================================

def calc_mass_mole_conversion(val: float, val_type: str, molar_mass: float, 
                               source_unit: str, target_unit: str) -> CalculationEntry:
    """Converts between mass and moles."""
    if val <= 0 or molar_mass <= 0:
        raise ValueError("Values and Molar Mass must be strictly greater than zero.")

    if val_type == 'mass_to_mole':
        # Convert input mass to grams
        mass_g = val * get_unit_factor(MASS_UNITS, source_unit, "mass")
        moles = mass_g / molar_mass
        res_val = moles / get_unit_factor(MOLE_UNITS, target_unit, "mole")
        expl = f"Moles = (Mass in g) / (Molar Mass in g/mol) = ({mass_g:.4f} g) / ({molar_mass} g/mol)"
        entry = CalculationEntry(
            operation="Mass -> Moles",
            inputs={"Mass": f"{val} {source_unit}", "Molar Mass": f"{molar_mass} g/mol"},
            result=f"{res_val:.10g} {target_unit}",
            explanation=expl
        )
    elif val_type == 'mole_to_mass':
        # Moles to mass
        moles = val * get_unit_factor(MOLE_UNITS, source_unit, "mole")
        mass_g = moles * molar_mass
        res_val = mass_g / get_unit_factor(MASS_UNITS, target_unit, "mass")
        expl = f"Mass = Moles * Molar Mass = ({moles:.4f} mol) * ({molar_mass} g/mol)"
        entry = CalculationEntry(
            operation="Moles -> Mass",
            inputs={"Moles": f"{val} {source_unit}", "Molar Mass": f"{molar_mass} g/mol"},
            result=f"{res_val:.10g} {target_unit}",
            explanation=expl
        )
    else:
        raise ValueError(
            "Invalid conversion type. Use 'mass_to_mole' or 'mole_to_mass'."
        )
    return entry


def solve_ideal_gas(target_var: str, p_tuple: Tuple[float, str], v_tuple: Tuple[float, str], 
                    n_tuple: Tuple[float, str], t_tuple: Tuple[float, str], 
                    out_unit: str) -> CalculationEntry:
    """
    Solves Ideal Gas Law (PV = nRT) for target_var ('P', 'V', 'n', 'T').
    """
    target = target_var.upper()
    if target not in ['P', 'V', 'N', 'T']:
        raise ValueError("Invalid target variable. Choose P, V, n, or T.")
    
    # Internal SI representations
    p_pa = p_tuple[0] * get_unit_factor(PRESSURE_UNITS, p_tuple[1], "pressure") if target != 'P' else None
    v_m3 = v_tuple[0] * get_unit_factor(VOLUME_UNITS, v_tuple[1], "volume") if target != 'V' else None
    n_mol = n_tuple[0] * get_unit_factor(MOLE_UNITS, n_tuple[1], "mole") if target != 'N' else None
    t_k = to_kelvin(t_tuple[0], t_tuple[1]) if target != 'T' else None
    
    # Physical validity checks
    if p_pa is not None and p_pa <= 0:
        raise ValueError("Pressure must be > 0.")
    if v_m3 is not None and v_m3 <= 0:
        raise ValueError("Volume must be > 0.")
    if n_mol is not None and n_mol <= 0:
        raise ValueError("Mole amount must be > 0.")
    if t_k is not None and t_k <= 0:
        raise ValueError("Temperature must be above absolute zero (> 0 K).")

    inputs_summary = {}
    if target != 'P': inputs_summary['Pressure'] = f"{p_tuple[0]} {p_tuple[1]}"
    if target != 'V': inputs_summary['Volume'] = f"{v_tuple[0]} {v_tuple[1]}"
    if target != 'N': inputs_summary['Moles'] = f"{n_tuple[0]} {n_tuple[1]}"
    if target != 'T': inputs_summary['Temperature'] = f"{t_tuple[0]} {t_tuple[1]}"

    if target == 'P':
        res_pa = (n_mol * R_SI * t_k) / v_m3
        res_final = res_pa / get_unit_factor(PRESSURE_UNITS, out_unit, "pressure")
        expl = f"P = (n * R * T) / V = ({n_mol:.4e} mol * {R_SI} J/(mol*K) * {t_k:.2f} K) / ({v_m3:.4e} m^3)"
        op = "Ideal Gas Solver (Find P)"
    elif target == 'V':
        res_m3 = (n_mol * R_SI * t_k) / p_pa
        res_final = res_m3 / get_unit_factor(VOLUME_UNITS, out_unit, "volume")
        expl = f"V = (n * R * T) / P = ({n_mol:.4e} mol * {R_SI} J/(mol*K) * {t_k:.2f} K) / ({p_pa:.4e} Pa)"
        op = "Ideal Gas Solver (Find V)"
    elif target == 'N':
        res_mol = (p_pa * v_m3) / (R_SI * t_k)
        res_final = res_mol / get_unit_factor(MOLE_UNITS, out_unit, "mole")
        expl = f"n = (P * V) / (R * T) = ({p_pa:.4e} Pa * {v_m3:.4e} m^3) / ({R_SI} J/(mol*K) * {t_k:.2f} K)"
        op = "Ideal Gas Solver (Find n)"
    elif target == 'T':
        res_k = (p_pa * v_m3) / (n_mol * R_SI)
        res_final = from_kelvin(res_k, out_unit)
        expl = f"T = (P * V) / (n * R) = ({p_pa:.4e} Pa * {v_m3:.4e} m^3) / ({n_mol:.4e} mol * {R_SI} J/(mol*K))"
        op = "Ideal Gas Solver (Find T)"
    return CalculationEntry(
        operation=op,
        inputs=inputs_summary,
        result=f"{res_final:.10g} {out_unit}",
        explanation=expl
    )


def calc_gas_volume_from_mass(mass_val: float, mass_unit: str, molar_mass: float,
                              p_val: float, p_unit: str, t_val: float, t_unit: str,
                              target_v_unit: str) -> CalculationEntry:
    """Calculates gas volume given mass, molar mass, temperature, and pressure."""
    if mass_val <= 0 or molar_mass <= 0:
        raise ValueError("Mass and Molar Mass must be strictly greater than zero.")

    # Step 1: Mass -> Moles
    mass_g = mass_val * get_unit_factor(MASS_UNITS, mass_unit, "mass")
    moles = mass_g / molar_mass

    # Step 2: Ideal Gas Law V = nRT / P
    calc_entry = solve_ideal_gas('V', (p_val, p_unit), (0, target_v_unit), (moles, 'mol'), (t_val, t_unit), target_v_unit)
    calc_entry.operation = "Gas Volume from Mass & State Parameters"
    calc_entry.inputs['Mass'] = f"{mass_val} {mass_unit}"
    calc_entry.inputs['Molar Mass'] = f"{molar_mass} g/mol"
    calc_entry.explanation = f"Calculated moles ({moles:.4f} mol) from mass, then solved " + calc_entry.explanation
    return calc_entry


def _require_positive(value: float, label: str) -> None:
    """Reject zero and negative physical quantities."""
    if value <= 0:
        raise ValueError(f"{label} must be greater than zero.")


def _volume_to_liters(value: float, unit: str) -> float:
    """Convert a volume to liters for solution calculations."""
    return value * get_unit_factor(VOLUME_UNITS, unit, "volume") * 1000.0


def calc_molarity(target_var: str, m_tuple: Tuple[float, str],
                  n_tuple: Tuple[float, str], v_tuple: Tuple[float, str],
                  out_unit: str) -> CalculationEntry:
    """Solve M = n/V for molarity, amount of substance, or volume."""
    target = target_var.upper()
    if target not in {'M', 'N', 'V'}:
        raise ValueError("Invalid molarity target. Choose M, n, or V.")

    molarity = (
        m_tuple[0] * get_unit_factor(MOLARITY_UNITS, m_tuple[1], "molarity")
        if target != 'M' else None
    )
    moles = (
        n_tuple[0] * get_unit_factor(MOLE_UNITS, n_tuple[1], "mole")
        if target != 'N' else None
    )
    volume_l = (
        _volume_to_liters(v_tuple[0], v_tuple[1])
        if target != 'V' else None
    )

    if molarity is not None:
        _require_positive(molarity, "Molarity")
    if moles is not None:
        _require_positive(moles, "Mole amount")
    if volume_l is not None:
        _require_positive(volume_l, "Volume")

    if target == 'M':
        result_base = moles / volume_l
        result = result_base / get_unit_factor(MOLARITY_UNITS, out_unit, "molarity")
        explanation = f"M = n / V = {moles:.6g} mol / {volume_l:.6g} L"
    elif target == 'N':
        result_base = molarity * volume_l
        result = result_base / get_unit_factor(MOLE_UNITS, out_unit, "mole")
        explanation = f"n = M * V = {molarity:.6g} mol/L * {volume_l:.6g} L"
    else:
        result_l = moles / molarity
        result = (result_l / 1000.0) / get_unit_factor(
            VOLUME_UNITS, out_unit, "volume"
        )
        explanation = f"V = n / M = {moles:.6g} mol / {molarity:.6g} mol/L"

    inputs = {
        key: f"{value[0]} {value[1]}"
        for key, value in zip(['M', 'n', 'V'], [m_tuple, n_tuple, v_tuple])
        if key != target
    }
    return CalculationEntry(
        operation=f"Molarity (Find {target})",
        inputs=inputs,
        result=f"{result:.10g} {out_unit}",
        explanation=explanation
    )


def calc_dilution(target_var: str, m1: Tuple[float, str], v1: Tuple[float, str],
                  m2: Tuple[float, str], v2: Tuple[float, str],
                  out_unit: str) -> CalculationEntry:
    """Solve M1V1 = M2V2 for any one dilution variable."""
    target = target_var.upper()
    if target not in {'M1', 'V1', 'M2', 'V2'}:
        raise ValueError("Invalid dilution target. Choose M1, V1, M2, or V2.")

    M1 = (
        m1[0] * get_unit_factor(MOLARITY_UNITS, m1[1], "molarity")
        if target != 'M1' else None
    )
    V1 = (
        _volume_to_liters(v1[0], v1[1])
        if target != 'V1' else None
    )
    M2 = (
        m2[0] * get_unit_factor(MOLARITY_UNITS, m2[1], "molarity")
        if target != 'M2' else None
    )
    V2 = (
        _volume_to_liters(v2[0], v2[1])
        if target != 'V2' else None
    )

    for value, label in [(M1, "M1"), (V1, "V1"), (M2, "M2"), (V2, "V2")]:
        if value is not None:
            _require_positive(value, label)

    if target == 'M1':
        result = (M2 * V2) / V1
        result /= get_unit_factor(MOLARITY_UNITS, out_unit, "molarity")
        explanation = f"M1 = (M2 * V2) / V1 = ({M2:.6g} * {V2:.6g}) / {V1:.6g}"
    elif target == 'V1':
        result_l = (M2 * V2) / M1
        result = (result_l / 1000.0) / get_unit_factor(
            VOLUME_UNITS, out_unit, "volume"
        )
        explanation = f"V1 = (M2 * V2) / M1 = ({M2:.6g} * {V2:.6g}) / {M1:.6g}"
    elif target == 'M2':
        result = (M1 * V1) / V2
        result /= get_unit_factor(MOLARITY_UNITS, out_unit, "molarity")
        explanation = f"M2 = (M1 * V1) / V2 = ({M1:.6g} * {V1:.6g}) / {V2:.6g}"
    else:
        result_l = (M1 * V1) / M2
        result = (result_l / 1000.0) / get_unit_factor(
            VOLUME_UNITS, out_unit, "volume"
        )
        explanation = f"V2 = (M1 * V1) / M2 = ({M1:.6g} * {V1:.6g}) / {M2:.6g}"

    inputs = {
        key: f"{value[0]} {value[1]}"
        for key, value in zip(['M1', 'V1', 'M2', 'V2'], [m1, v1, m2, v2])
        if key != target
    }
    return CalculationEntry(
        operation=f"Dilution (Find {target})",
        inputs=inputs,
        result=f"{result:.10g} {out_unit}",
        explanation=explanation
    )


def _volume_to_milliliters(value: float, unit: str) -> float:
    """Convert a volume to milliliters for density calculations."""
    return value * get_unit_factor(VOLUME_UNITS, unit, "volume") * 1e6


def calc_liquid_density(target_var: str, mass: Tuple[float, str],
                        volume: Tuple[float, str], density: Tuple[float, str],
                        out_unit: str) -> CalculationEntry:
    """Solve d = m/V for liquid density, mass, or volume."""
    target = target_var.upper()
    if target not in {'D', 'M', 'V'}:
        raise ValueError("Invalid density target. Choose d, m, or V.")

    mass_g = (
        mass[0] * get_unit_factor(MASS_UNITS, mass[1], "mass")
        if target != 'M' else None
    )
    volume_ml = (
        _volume_to_milliliters(volume[0], volume[1])
        if target != 'V' else None
    )
    density_gml = (
        density[0] * get_unit_factor(DENSITY_UNITS, density[1], "density")
        if target != 'D' else None
    )

    for value, label in [
        (mass_g, "Mass"), (volume_ml, "Volume"), (density_gml, "Density")
    ]:
        if value is not None:
            _require_positive(value, label)

    if target == 'D':
        result = (mass_g / volume_ml) / get_unit_factor(
            DENSITY_UNITS, out_unit, "density"
        )
        explanation = f"d = m / V = {mass_g:.6g} g / {volume_ml:.6g} mL"
    elif target == 'M':
        result = (density_gml * volume_ml) / get_unit_factor(
            MASS_UNITS, out_unit, "mass"
        )
        explanation = f"m = d * V = {density_gml:.6g} g/mL * {volume_ml:.6g} mL"
    else:
        result_ml = mass_g / density_gml
        result = (result_ml / 1e6) / get_unit_factor(
            VOLUME_UNITS, out_unit, "volume"
        )
        explanation = f"V = m / d = {mass_g:.6g} g / {density_gml:.6g} g/mL"

    inputs = {
        key: f"{value[0]} {value[1]}"
        for key, value in zip(['m', 'V', 'd'], [mass, volume, density])
        if key != target
    }
    return CalculationEntry(
        operation=f"Liquid Density (Find {target})",
        inputs=inputs,
        result=f"{result:.10g} {out_unit}",
        explanation=explanation
    )


# ==============================================================================
# BATCH FILE PROCESSING ENGINE
# ==============================================================================

def process_batch_file(filepath: str) -> None:
    """
    Reads a CSV file containing batch jobs.
    CSV Format expected:
    job_type, arg1, arg2, arg3, arg4, arg5, arg6, arg7, target_unit

    Supported job types:
      mass_to_mole, mole_to_mass, ideal_gas, molarity, dilution, density
    """
    if not os.path.exists(filepath):
        print(f"Error: Batch file '{filepath}' not found.")
        return

    print(f"\n--- PROCESSING BATCH FILE: {filepath} ---")
    with open(filepath, mode='r', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = next(reader, None) # Skip header line
        
        for idx, row in enumerate(reader, start=1):
            if not row or row[0].startswith('#'):
                continue
            try:
                job_type = row[0].strip().lower()
                if job_type == 'mass_to_mole':
                    res = calc_mass_mole_conversion(
                        float(row[1]), 'mass_to_mole', calculate_molar_mass(row[2]), row[3], row[4]
                    )
                elif job_type == 'mole_to_mass':
                    res = calc_mass_mole_conversion(
                        float(row[1]), 'mole_to_mass', calculate_molar_mass(row[2]), row[3], row[4]
                    )
                elif job_type == 'ideal_gas':
                    target = row[1].strip()
                    res = solve_ideal_gas(
                        target,
                        parse_environment_input(row[2], row[3], "pressure"), # P
                        (float(row[4]), row[5]), # V
                        (float(row[6]), row[7]), # n
                        parse_environment_input(row[8], row[9], "temperature"), # T
                        row[10].strip()          # output unit
                    )
                elif job_type == 'molarity':
                    res = calc_molarity(
                        row[1].strip(),
                        (float(row[2]), row[3]),  # M
                        (float(row[4]), row[5]),  # n
                        (float(row[6]), row[7]),  # V
                        row[8].strip()            # output unit
                    )
                elif job_type == 'dilution':
                    res = calc_dilution(
                        row[1].strip(),
                        (float(row[2]), row[3]),  # M1
                        (float(row[4]), row[5]),  # V1
                        (float(row[6]), row[7]),  # M2
                        (float(row[8]), row[9]),  # V2
                        row[10].strip()           # output unit
                    )
                elif job_type in {'density', 'liquid_density'}:
                    res = calc_liquid_density(
                        row[1].strip(),
                        (float(row[2]), row[3]),  # mass
                        (float(row[4]), row[5]),  # volume
                        (float(row[6]), row[7]),  # density
                        row[8].strip()            # output unit
                    )
                else:
                    print(f"Row {idx}: Unknown job_type '{job_type}'")
                    continue

                HISTORY.append(res)
                print(f"Job #{idx} Success -> {res.operation}: Result = {res.result}")
            except Exception as e:
                print(f"Row {idx} Failed: {e}")
    print("--- BATCH PROCESSING COMPLETE ---\n")


def generate_sample_csv(filename: str = "batch_sample.csv") -> None:
    """Generates a sample batch CSV file for demonstration."""
    sample_data = [
        ["# JobType", "Param1", "Param2", "Param3", "Param4", "Param5", "Param6", "Param7", "Param8", "Param9", "TargetUnit"],
        ["mass_to_mole", "18.015", "18.015", "g", "mol"],
        ["mole_to_mass", "2.5", "32.0", "mol", "g"],
        ["ideal_gas", "V", "1.0", "atm", "0", "L", "2.0", "mol", "25", "C", "L"]
    ]
    with open(filename, mode='w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerows(sample_data)
    print(f"Created sample batch template at '{filename}'")


# ==============================================================================
# INTERACTIVE CLI INTERFACE
# ==============================================================================

def prompt_float(prompt_text: str) -> float:
    while True:
        try:
            return float(input(prompt_text).strip())
        except ValueError:
            print("Invalid input. Please enter a numerical value.")

def prompt_string(prompt_text: str) -> str:
    return input(prompt_text).strip()

def prompt_molar_mass() -> float:
    """Prompt for a numeric molar mass or a chemical formula."""
    while True:
        raw_value = prompt_string(
            "Enter Molar Mass in g/mol (number or formula such as H2O, CuSO4): "
        )
        try:
            mass = calculate_molar_mass(raw_value)
            print(f"  -> Using Molar Mass: {mass:.4f} g/mol")
            return mass
        except ValueError as err:
            print(f"Invalid molar mass: {err}")

def prompt_environment(env_type: str) -> Tuple[float, str]:
    """
    Prompt for temperature or pressure, supporting common presets.

    Temperature presets:
      STP = 0 C, RT/SATP = 25 C
    Pressure presets:
      STP/RT/SATP = 1 atm
    """
    if env_type not in {"Temperature", "Pressure"}:
        raise ValueError("Environment type must be Temperature or Pressure.")

    while True:
        if env_type == "Temperature":
            raw_value = prompt_string(
                "Temperature value or preset (STP=0 C, RT/SATP=25 C): "
            ).upper()
            if raw_value == "STP":
                return 0.0, "C"
            if raw_value in {"RT", "SATP"}:
                return 25.0, "C"
            try:
                return float(raw_value), prompt_string(
                    "Temperature unit (C, K, F): "
                )
            except ValueError:
                print("Invalid temperature. Enter a number, STP, RT, or SATP.")
        else:
            raw_value = prompt_string(
                "Pressure value or preset (STP/RT/SATP=1 atm): "
            ).upper()
            if raw_value in {"STP", "RT", "SATP"}:
                return 1.0, "atm"
            try:
                return float(raw_value), prompt_string(
                    "Pressure unit (atm, Pa, kPa, Torr): "
                )
            except ValueError:
                print("Invalid pressure. Enter a number, STP, RT, or SATP.")

def show_unit_reference_table():
    print("\n" + "="*55)
    print("            SUPPORTED UNITS REFERENCE TABLE")
    print("="*55)
    print("Mass:        g, kg, mg, ug, lb, oz")
    print("Volume:      L, mL, m3, cm3, cc, gal")
    print("Moles:       mol, mmol, kmol, umol")
    print("Pressure:    atm, Pa, kPa, MPa, bar, mbar, Torr, mmHg, psi")
    print("Temperature: C, K, F")
    print("="*55 + "\n")

def display_entry(entry: CalculationEntry):
    print("\n" + "-"*50)
    print(f" OPERATION:   {entry.operation}")
    print(" INPUTS:")
    for k, v in entry.inputs.items():
        print(f"   * {k}: {v}")
    print(f" RESULT:      {entry.result}")
    print(f" FORMULA/EXPLANATION:")
    print(f"   {entry.explanation}")
    print("-" * 50 + "\n")

def main_menu():
    while True:
        print("==================================================")
        print("   MASS, VOLUME, AND GAS MOLE CALCULATOR")
        print("==================================================")
        print("1. Convert Mass <-> Moles (formula molar masses supported)")
        print("2. Solve Ideal Gas Law (PV = nRT; STP/RT presets supported)")
        print("3. Calculate Gas Volume from Mass & State (P, T; presets supported)")
        print("4. Quick Unit Converter")
        print("5. Run Batch Processing from CSV File")
        print("6. Generate Sample Batch CSV File")
        print("7. View Session History Log")
        print("8. Display Supported Units Table")
        print("9. Exit")
        print("10. Solve Molarity (M = n / V)")
        print("11. Solve Dilution (M1V1 = M2V2)")
        print("12. Solve Liquid Density (d = m / V)")
        print("==================================================")
        
        choice = prompt_string("Select an option (1-12): ")
        
        try:
            if choice == '1':
                sub = prompt_string("Choose: (1) Mass to Moles, (2) Moles to Mass: ")
                if sub not in {'1', '2'}:
                    raise ValueError("Choose either 1 for Mass to Moles or 2 for Moles to Mass.")
                val = prompt_float("Enter input value: ")
                src_u = prompt_string("Enter input unit (e.g., g, kg, mol): ")
                mm = prompt_molar_mass()
                tgt_u = prompt_string("Enter desired output unit (e.g., mol, mmol, g): ")
                
                job = 'mass_to_mole' if sub == '1' else 'mole_to_mass'
                entry = calc_mass_mole_conversion(val, job, mm, src_u, tgt_u)
                HISTORY.append(entry)
                display_entry(entry)

            elif choice == '2':
                target = prompt_string("What variable do you want to solve for? (P, V, n, T): ").upper()
                if target not in ['P', 'V', 'N', 'T']:
                    print("Invalid target choice.")
                    continue

                p_tuple = prompt_environment("Pressure") if target != 'P' else (0.0, '')
                v_tuple = (prompt_float("Volume value: "), prompt_string("Volume unit (L, mL, m3): ")) if target != 'V' else (0.0, '')
                n_tuple = (prompt_float("Moles value: "), prompt_string("Mole unit (mol, mmol): ")) if target != 'N' else (0.0, '')
                t_tuple = prompt_environment("Temperature") if target != 'T' else (0.0, '')
                
                out_u = prompt_string(f"Enter target unit for {target}: ")
                entry = solve_ideal_gas(target, p_tuple, v_tuple, n_tuple, t_tuple, out_u)
                HISTORY.append(entry)
                display_entry(entry)

            elif choice == '3':
                m_val = prompt_float("Enter Mass value: ")
                m_unit = prompt_string("Enter Mass unit (g, kg, lb): ")
                mm = prompt_molar_mass()
                p_val, p_unit = prompt_environment("Pressure")
                t_val, t_unit = prompt_environment("Temperature")
                v_out_unit = prompt_string("Enter desired Volume unit (L, mL, m3): ")
                
                entry = calc_gas_volume_from_mass(m_val, m_unit, mm, p_val, p_unit, t_val, t_unit, v_out_unit)
                HISTORY.append(entry)
                display_entry(entry)

            elif choice == '4':
                cat = prompt_string("Category ((M)ass, (V)olume, (P)ressure, (MOL)e): ").upper()
                val = prompt_float("Enter value: ")
                src = prompt_string("From unit: ")
                tgt = prompt_string("To unit: ")
                
                if cat.startswith('M') and not cat.startswith('MOL'):
                    res = convert_mass(val, src, tgt)
                elif cat.startswith('V'):
                    res = convert_volume(val, src, tgt)
                elif cat.startswith('P'):
                    res = convert_pressure(val, src, tgt)
                elif cat.startswith('MOL'):
                    res = convert_mole(val, src, tgt)
                else:
                    print("Unknown category.")
                    continue
                print(f"\nConversions Result: {val} {src} = {res:.10g} {tgt}\n")

            elif choice == '5':
                fpath = prompt_string("Enter CSV file path: ")
                process_batch_file(fpath)

            elif choice == '6':
                generate_sample_csv()

            elif choice == '7':
                if not HISTORY:
                    print("\nNo calculations recorded in this session yet.\n")
                else:
                    print(f"\n=== SESSION HISTORY ({len(HISTORY)} entries) ===")
                    for idx, item in enumerate(HISTORY, start=1):
                        print(f"[{idx}] {item.operation} -> Result: {item.result}")
                    print("==========================================\n")

            elif choice == '8':
                show_unit_reference_table()

            elif choice == '10':
                target = prompt_string("Solve molarity for (M, n, V): ").upper()
                if target not in {'M', 'N', 'V'}:
                    raise ValueError("Choose M, n, or V.")

                m_tuple = (
                    prompt_float("Molarity value: "),
                    prompt_string("Molarity unit (M, mol/L, mM): ")
                ) if target != 'M' else (0.0, '')

                n_tuple = (0.0, '')
                if target != 'N':
                    n_mode = prompt_string(
                        "Mole input: (1) amount in moles, (2) mass + formula: "
                    )
                    if n_mode == '1':
                        n_tuple = (
                            prompt_float("Moles value: "),
                            prompt_string("Mole unit (mol, mmol): ")
                        )
                    elif n_mode == '2':
                        mass_value = prompt_float("Mass value: ")
                        mass_unit = prompt_string("Mass unit (g, kg, etc.): ")
                        molar_mass = prompt_molar_mass()
                        converted = calc_mass_mole_conversion(
                            mass_value, 'mass_to_mole', molar_mass,
                            mass_unit, 'mol'
                        )
                        print(f"  -> Converted to {converted.result}")
                        n_tuple = (float(converted.result.split()[0]), 'mol')
                    else:
                        raise ValueError("Choose mole input option 1 or 2.")

                v_tuple = (
                    prompt_float("Volume value: "),
                    prompt_string("Volume unit (L, mL, m3): ")
                ) if target != 'V' else (0.0, '')

                out_unit = prompt_string(f"Target unit for {target}: ")
                entry = calc_molarity(
                    target, m_tuple, n_tuple, v_tuple, out_unit
                )
                HISTORY.append(entry)
                display_entry(entry)

            elif choice == '11':
                target = prompt_string(
                    "Solve dilution for (M1, V1, M2, V2): "
                ).upper()
                if target not in {'M1', 'V1', 'M2', 'V2'}:
                    raise ValueError("Choose M1, V1, M2, or V2.")

                m1 = (
                    prompt_float("M1 value: "),
                    prompt_string("M1 unit (M, mol/L, mM): ")
                ) if target != 'M1' else (0.0, '')
                v1 = (
                    prompt_float("V1 value: "),
                    prompt_string("V1 unit (L, mL, m3): ")
                ) if target != 'V1' else (0.0, '')
                m2 = (
                    prompt_float("M2 value: "),
                    prompt_string("M2 unit (M, mol/L, mM): ")
                ) if target != 'M2' else (0.0, '')
                v2 = (
                    prompt_float("V2 value: "),
                    prompt_string("V2 unit (L, mL, m3): ")
                ) if target != 'V2' else (0.0, '')
                out_unit = prompt_string(f"Target unit for {target}: ")
                entry = calc_dilution(
                    target, m1, v1, m2, v2, out_unit
                )
                HISTORY.append(entry)
                display_entry(entry)

            elif choice == '12':
                target = prompt_string(
                    "Solve density for (d, m, V): "
                ).upper()
                if target not in {'D', 'M', 'V'}:
                    raise ValueError("Choose d, m, or V.")

                mass = (
                    prompt_float("Mass value: "),
                    prompt_string("Mass unit (g, kg, etc.): ")
                ) if target != 'M' else (0.0, '')
                volume = (
                    prompt_float("Volume value: "),
                    prompt_string("Volume unit (mL, L, m3): ")
                ) if target != 'V' else (0.0, '')
                density = (
                    prompt_float("Density value: "),
                    prompt_string("Density unit (g/mL, kg/L): ")
                ) if target != 'D' else (0.0, '')
                out_unit = prompt_string(f"Target unit for {target}: ")
                entry = calc_liquid_density(
                    target, mass, volume, density, out_unit
                )
                HISTORY.append(entry)
                display_entry(entry)

            elif choice == '9':
                print("Exiting calculator. Goodbye!")
                sys.exit(0)

            else:
                print("Invalid option. Please choose 1-12.")

        except Exception as err:
            print(f"\n[ERROR]: {err}\n")


if __name__ == '__main__':
    # Parse potential command line batch execution flag
    if len(sys.argv) > 1:
        if sys.argv[1] in ['--batch', '-b'] and len(sys.argv) > 2:
            process_batch_file(sys.argv[2])
        elif sys.argv[1] in ['--help', '-h']:
            print("Usage: python3 script.py [--batch filepath.csv]")
        else:
            print(f"Unknown argument. Running interactive menu...")
            main_menu()
    else:
        main_menu()
