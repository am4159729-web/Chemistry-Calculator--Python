#!/usr/bin/env python3
"""Mass, mole, gas-law and solution calculator with an interactive menu and a CSV batch mode."""
import argparse
import csv
import math
import re
from dataclasses import dataclass

R = 8.31446261815324  # J/(mol*K)

# Standard atomic weights, g/mol.
ATOMIC_MASSES = {
    "H": 1.008, "He": 4.0026, "Li": 6.94, "Be": 9.0122, "B": 10.81,
    "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.180,
    "Na": 22.989, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974,
    "S": 32.06, "Cl": 35.45, "Ar": 39.95, "K": 39.098, "Ca": 40.078,
    "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996, "Mn": 54.938,
    "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.38,
    "Ga": 69.723, "Ge": 72.630, "As": 74.922, "Se": 78.971, "Br": 79.904,
    "Kr": 83.798, "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224,
    "Nb": 92.906, "Mo": 95.95, "Ru": 101.07, "Rh": 102.91, "Pd": 106.42,
    "Ag": 107.87, "Cd": 112.41, "In": 114.82, "Sn": 118.71, "Sb": 121.76,
    "Te": 127.60, "I": 126.90, "Xe": 131.29, "Cs": 132.91, "Ba": 137.33,
    "W": 183.84, "Pt": 195.08, "Au": 196.97, "Hg": 200.59, "Pb": 207.2,
    "Bi": 208.98, "U": 238.03,
}

# Unit -> factor to the base unit: g, m^3, mol, Pa, mol/L, g/mL.
UNITS = {
    "mass": {
        "g": 1.0, "gram": 1.0, "grams": 1.0,
        "kg": 1000.0, "kilogram": 1000.0, "kilograms": 1000.0,
        "mg": 0.001, "milligram": 0.001, "milligrams": 0.001,
        "ug": 1e-6, "microgram": 1e-6, "micrograms": 1e-6,
        "lb": 453.59237, "lbs": 453.59237, "pound": 453.59237, "pounds": 453.59237,
        "oz": 28.349523125, "ounce": 28.349523125, "ounces": 28.349523125,
    },
    "volume": {
        "m3": 1.0, "m^3": 1.0, "cubic_meter": 1.0, "cubic_meters": 1.0,
        "l": 0.001, "liter": 0.001, "liters": 0.001, "litre": 0.001, "litres": 0.001,
        "ml": 1e-6, "milliliter": 1e-6, "milliliters": 1e-6, "cc": 1e-6, "cm3": 1e-6, "cm^3": 1e-6,
        "gal": 0.003785411784, "gallon": 0.003785411784, "gallons": 0.003785411784,
    },
    "mole": {
        "mol": 1.0, "mole": 1.0, "moles": 1.0,
        "mmol": 0.001, "millimole": 0.001, "millimoles": 0.001,
        "kmol": 1000.0, "kilomole": 1000.0, "kilomoles": 1000.0,
        "umol": 1e-6, "micromole": 1e-6, "micromoles": 1e-6,
    },
    "pressure": {
        "pa": 1.0, "pascal": 1.0, "pascals": 1.0,
        "kpa": 1000.0, "kilopascal": 1000.0, "kilopascals": 1000.0,
        "mpa": 1e6, "megapascal": 1e6, "megapascals": 1e6,
        "atm": 101325.0, "atmosphere": 101325.0, "atmospheres": 101325.0,
        "torr": 101325.0 / 760.0, "torrs": 101325.0 / 760.0, "mmhg": 101325.0 / 760.0,
        "bar": 100000.0, "bars": 100000.0,
        "mbar": 100.0, "millibar": 100.0,
        "psi": 6894.757293168,
    },
    "molarity": {
        "m": 1.0, "molar": 1.0, "mol/l": 1.0, "mol/dm3": 1.0,
        "mm": 0.001, "millimolar": 0.001, "mmol/l": 0.001,
    },
    "density": {
        "g/ml": 1.0, "g/cm3": 1.0, "g/cc": 1.0, "kg/l": 1.0,
        "kg/m3": 0.001, "g/l": 0.001,
    },
}

# STP, RT and SATP all mean 1 atm; the value is the preset temperature in deg C.
PRESET_TEMPS_C = {"STP": 0.0, "RT": 25.0, "SATP": 25.0}

Quantity = tuple[float, str]  # (value, unit)
UNKNOWN: Quantity = (0.0, "")  # placeholder for the variable being solved for

MENU = """\
1.  Mass <-> moles
2.  Ideal gas law (PV = nRT)
3.  Gas volume from mass
4.  Unit converter
5.  Run batch CSV file
6.  Write sample batch CSV
7.  Session history
8.  Supported units
9.  Exit
10. Molarity (M = n/V)
11. Dilution (M1V1 = M2V2)
12. Liquid density (d = m/V)"""

UNIT_HELP = """\
Mass:        g, kg, mg, ug, lb, oz
Volume:      L, mL, m3, cm3, cc, gal
Moles:       mol, mmol, kmol, umol
Pressure:    atm, Pa, kPa, MPa, bar, mbar, Torr, mmHg, psi
Temperature: C, K, F
Molarity:    M, mM, mol/L, mmol/L
Density:     g/mL, g/cm3, g/cc, kg/L, kg/m3, g/L"""


@dataclass
class CalculationEntry:
    operation: str
    inputs: dict[str, str]
    result: str
    explanation: str

    def __str__(self) -> str:
        lines = [self.operation]
        lines += [f"  {name}: {value}" for name, value in self.inputs.items()]
        lines += [f"  {self.explanation}", f"  Result: {self.result}"]
        return "\n".join(lines)


HISTORY: list[CalculationEntry] = []


def unit_factor(kind: str, unit: str) -> float:
    table = UNITS[kind]
    try:
        return table[unit.strip().lower()]
    except KeyError:
        supported = ", ".join(sorted(table))
        raise ValueError(
            f"Unsupported {kind} unit '{unit}'. Supported units: {supported}"
        ) from None


def convert(kind: str, value: float, src: str, dst: str) -> float:
    return value * unit_factor(kind, src) / unit_factor(kind, dst)


def to_kelvin(value: float, unit: str) -> float:
    u = unit.strip().lower()
    if u in ("k", "kelvin"):
        return value
    if u in ("c", "celsius"):
        return value + 273.15
    if u in ("f", "fahrenheit"):
        return (value - 32.0) * (5.0 / 9.0) + 273.15
    raise ValueError(f"Unsupported temperature unit '{unit}'. Supported units: K, C, F")


def from_kelvin(kelvin: float, unit: str) -> float:
    u = unit.strip().lower()
    if u in ("k", "kelvin"):
        return kelvin
    if u in ("c", "celsius"):
        return kelvin - 273.15
    if u in ("f", "fahrenheit"):
        return (kelvin - 273.15) * (9.0 / 5.0) + 32.0
    raise ValueError(f"Unsupported temperature unit '{unit}'. Supported units: K, C, F")


def require_positive(**values) -> None:
    for name, value in values.items():
        if value is not None and value <= 0:
            raise ValueError(f"{name.replace('_', ' ').capitalize()} must be greater than zero.")


def calculate_molar_mass(text: str) -> float:
    """g/mol from a number or a flat formula such as CuSO4 (no parentheses or hydrates)."""
    text = text.strip()
    try:
        mass = float(text)
    except ValueError:
        pass
    else:
        if not math.isfinite(mass) or mass <= 0:
            raise ValueError("Molar mass must be a finite value greater than zero.")
        return mass

    if not re.fullmatch(r"(?:[A-Z][a-z]*\d*)+", text):
        raise ValueError(
            f"Could not parse '{text}'. Enter a positive number or a formula such as H2O or CuSO4."
        )
    total = 0.0
    for element, digits in re.findall(r"([A-Z][a-z]*)(\d*)", text):
        if element not in ATOMIC_MASSES:
            raise ValueError(f"Unknown element '{element}'. Check the chemical symbol.")
        count = int(digits) if digits else 1
        if count == 0:
            raise ValueError("Element counts in a formula must be greater than zero.")
        total += ATOMIC_MASSES[element] * count
    return total


def parse_environment_input(value: str, unit: str, env_type: str) -> Quantity:
    preset = value.strip().upper()
    if preset in PRESET_TEMPS_C:
        return (PRESET_TEMPS_C[preset], "C") if env_type == "temperature" else (1.0, "atm")
    try:
        return float(value), unit
    except ValueError:
        raise ValueError(f"Invalid {env_type} '{value}'. Use a number, STP, RT, or SATP.") from None


def calc_mass_mole_conversion(value: float, mode: str, molar_mass: float,
                              src_unit: str, dst_unit: str) -> CalculationEntry:
    require_positive(value=value, molar_mass=molar_mass)

    if mode == "mass_to_mole":
        mass_g = value * unit_factor("mass", src_unit)
        result = mass_g / molar_mass / unit_factor("mole", dst_unit)
        operation, label = "Mass -> Moles", "Mass"
        explanation = f"n = m / M = {mass_g:.6g} g / {molar_mass} g/mol"
    elif mode == "mole_to_mass":
        moles = value * unit_factor("mole", src_unit)
        result = moles * molar_mass / unit_factor("mass", dst_unit)
        operation, label = "Moles -> Mass", "Moles"
        explanation = f"m = n * M = {moles:.6g} mol * {molar_mass} g/mol"
    else:
        raise ValueError("Mode must be 'mass_to_mole' or 'mole_to_mass'.")

    inputs = {label: f"{value:.10g} {src_unit}", "Molar Mass": f"{molar_mass} g/mol"}
    return CalculationEntry(operation, inputs, f"{result:.10g} {dst_unit}", explanation)


def solve_ideal_gas(target_var: str, pressure: Quantity, volume: Quantity, moles: Quantity,
                    temperature: Quantity, out_unit: str) -> CalculationEntry:
    target = target_var.upper()
    if target not in ("P", "V", "N", "T"):
        raise ValueError("Invalid target variable. Choose P, V, n, or T.")

    p = None if target == "P" else pressure[0] * unit_factor("pressure", pressure[1])
    v = None if target == "V" else volume[0] * unit_factor("volume", volume[1])
    n = None if target == "N" else moles[0] * unit_factor("mole", moles[1])
    t = None if target == "T" else to_kelvin(*temperature)
    require_positive(pressure=p, volume=v, moles=n)
    if t is not None and t <= 0:
        raise ValueError("Temperature must be above absolute zero (> 0 K).")

    if target == "P":
        result = n * R * t / v / unit_factor("pressure", out_unit)
        explanation = f"P = nRT/V = ({n:.4e} mol * {R} J/(mol*K) * {t:.2f} K) / ({v:.4e} m^3)"
    elif target == "V":
        result = n * R * t / p / unit_factor("volume", out_unit)
        explanation = f"V = nRT/P = ({n:.4e} mol * {R} J/(mol*K) * {t:.2f} K) / ({p:.4e} Pa)"
    elif target == "N":
        result = p * v / (R * t) / unit_factor("mole", out_unit)
        explanation = f"n = PV/RT = ({p:.4e} Pa * {v:.4e} m^3) / ({R} J/(mol*K) * {t:.2f} K)"
    else:
        result = from_kelvin(p * v / (n * R), out_unit)
        explanation = f"T = PV/nR = ({p:.4e} Pa * {v:.4e} m^3) / ({n:.4e} mol * {R} J/(mol*K))"

    given = {"P": pressure, "V": volume, "n": moles, "T": temperature}
    inputs = {k: f"{val:.10g} {unit}" for k, (val, unit) in given.items() if k.upper() != target}
    return CalculationEntry(f"Ideal Gas Law (Find {target})", inputs,
                            f"{result:.10g} {out_unit}", explanation)


def calc_gas_volume_from_mass(mass: Quantity, molar_mass: float, pressure: Quantity,
                              temperature: Quantity, out_unit: str) -> CalculationEntry:
    require_positive(mass=mass[0], molar_mass=molar_mass)
    moles = mass[0] * unit_factor("mass", mass[1]) / molar_mass

    entry = solve_ideal_gas("V", pressure, UNKNOWN, (moles, "mol"), temperature, out_unit)
    entry.operation = "Gas Volume from Mass"
    entry.inputs = {"Mass": f"{mass[0]:.10g} {mass[1]}", "Molar Mass": f"{molar_mass} g/mol",
                    **entry.inputs}
    entry.explanation = f"n = {moles:.6g} mol from mass; {entry.explanation}"
    return entry


def calc_molarity(target_var: str, molarity: Quantity, moles: Quantity, volume: Quantity,
                  out_unit: str) -> CalculationEntry:
    target = target_var.upper()
    if target not in ("M", "N", "V"):
        raise ValueError("Invalid molarity target. Choose M, n, or V.")

    m = None if target == "M" else molarity[0] * unit_factor("molarity", molarity[1])
    n = None if target == "N" else moles[0] * unit_factor("mole", moles[1])
    liters = None if target == "V" else convert("volume", volume[0], volume[1], "L")
    require_positive(molarity=m, moles=n, volume=liters)

    if target == "M":
        result = n / liters / unit_factor("molarity", out_unit)
        explanation = f"M = n / V = {n:.6g} mol / {liters:.6g} L"
    elif target == "N":
        result = m * liters / unit_factor("mole", out_unit)
        explanation = f"n = M * V = {m:.6g} mol/L * {liters:.6g} L"
    else:
        result = convert("volume", n / m, "L", out_unit)
        explanation = f"V = n / M = {n:.6g} mol / {m:.6g} mol/L"

    given = {"M": molarity, "n": moles, "V": volume}
    inputs = {k: f"{val:.10g} {unit}" for k, (val, unit) in given.items() if k.upper() != target}
    return CalculationEntry(f"Molarity (Find {target})", inputs,
                            f"{result:.10g} {out_unit}", explanation)


def calc_dilution(target_var: str, m1: Quantity, v1: Quantity, m2: Quantity, v2: Quantity,
                  out_unit: str) -> CalculationEntry:
    target = target_var.upper()
    if target not in ("M1", "V1", "M2", "V2"):
        raise ValueError("Invalid dilution target. Choose M1, V1, M2, or V2.")

    conc1 = None if target == "M1" else m1[0] * unit_factor("molarity", m1[1])
    vol1 = None if target == "V1" else convert("volume", v1[0], v1[1], "L")
    conc2 = None if target == "M2" else m2[0] * unit_factor("molarity", m2[1])
    vol2 = None if target == "V2" else convert("volume", v2[0], v2[1], "L")
    require_positive(m1=conc1, v1=vol1, m2=conc2, v2=vol2)

    if target == "M1":
        result = conc2 * vol2 / vol1 / unit_factor("molarity", out_unit)
        explanation = f"M1 = (M2 * V2) / V1 = ({conc2:.6g} * {vol2:.6g}) / {vol1:.6g}"
    elif target == "V1":
        result = convert("volume", conc2 * vol2 / conc1, "L", out_unit)
        explanation = f"V1 = (M2 * V2) / M1 = ({conc2:.6g} * {vol2:.6g}) / {conc1:.6g}"
    elif target == "M2":
        result = conc1 * vol1 / vol2 / unit_factor("molarity", out_unit)
        explanation = f"M2 = (M1 * V1) / V2 = ({conc1:.6g} * {vol1:.6g}) / {vol2:.6g}"
    else:
        result = convert("volume", conc1 * vol1 / conc2, "L", out_unit)
        explanation = f"V2 = (M1 * V1) / M2 = ({conc1:.6g} * {vol1:.6g}) / {conc2:.6g}"

    given = {"M1": m1, "V1": v1, "M2": m2, "V2": v2}
    inputs = {k: f"{val:.10g} {unit}" for k, (val, unit) in given.items() if k != target}
    return CalculationEntry(f"Dilution (Find {target})", inputs,
                            f"{result:.10g} {out_unit}", explanation)


def calc_liquid_density(target_var: str, mass: Quantity, volume: Quantity, density: Quantity,
                        out_unit: str) -> CalculationEntry:
    target = target_var.upper()
    if target not in ("D", "M", "V"):
        raise ValueError("Invalid density target. Choose d, m, or V.")

    m = None if target == "M" else mass[0] * unit_factor("mass", mass[1])
    ml = None if target == "V" else convert("volume", volume[0], volume[1], "mL")
    d = None if target == "D" else density[0] * unit_factor("density", density[1])
    require_positive(mass=m, volume=ml, density=d)

    if target == "D":
        result = m / ml / unit_factor("density", out_unit)
        explanation = f"d = m / V = {m:.6g} g / {ml:.6g} mL"
    elif target == "M":
        result = d * ml / unit_factor("mass", out_unit)
        explanation = f"m = d * V = {d:.6g} g/mL * {ml:.6g} mL"
    else:
        result = convert("volume", m / d, "mL", out_unit)
        explanation = f"V = m / d = {m:.6g} g / {d:.6g} g/mL"

    given = {"m": mass, "V": volume, "d": density}
    inputs = {k: f"{val:.10g} {unit}" for k, (val, unit) in given.items() if k.upper() != target}
    return CalculationEntry(f"Liquid Density (Find {target})", inputs,
                            f"{result:.10g} {out_unit}", explanation)


def row_quantity(row: list[str], i: int) -> Quantity:
    return float(row[i]), row[i + 1]


def process_batch_file(filepath: str) -> None:
    """Run every row of a CSV file. The first row is a header; rows starting with '#' are skipped.

    Columns by job type. The unknown variable's value/unit columns are placeholders
    but the value must still be numeric.
      mass_to_mole, mole_to_mass: job, value, molar mass or formula, from unit, to unit
      ideal_gas: job, target, P, unit, V, unit, n, unit, T, unit, output unit
                 (P and T also accept STP, RT or SATP; the unit column is ignored then)
      molarity:  job, target, M, unit, n, unit, V, unit, output unit
      dilution:  job, target, M1, unit, V1, unit, M2, unit, V2, unit, output unit
      density:   job, target, m, unit, V, unit, d, unit, output unit
    """
    with open(filepath, newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))

    for i, row in enumerate(rows[1:], start=1):
        row = [cell.strip() for cell in row]
        if not row or row[0].startswith("#"):
            continue
        job = row[0].lower()
        try:
            if job in ("mass_to_mole", "mole_to_mass"):
                entry = calc_mass_mole_conversion(
                    float(row[1]), job, calculate_molar_mass(row[2]), row[3], row[4])
            elif job == "ideal_gas":
                entry = solve_ideal_gas(
                    row[1],
                    parse_environment_input(row[2], row[3], "pressure"),
                    row_quantity(row, 4),
                    row_quantity(row, 6),
                    parse_environment_input(row[8], row[9], "temperature"),
                    row[10])
            elif job == "molarity":
                entry = calc_molarity(
                    row[1], row_quantity(row, 2), row_quantity(row, 4), row_quantity(row, 6),
                    row[8])
            elif job == "dilution":
                entry = calc_dilution(
                    row[1], row_quantity(row, 2), row_quantity(row, 4), row_quantity(row, 6),
                    row_quantity(row, 8), row[10])
            elif job in ("density", "liquid_density"):
                entry = calc_liquid_density(
                    row[1], row_quantity(row, 2), row_quantity(row, 4), row_quantity(row, 6),
                    row[8])
            else:
                raise ValueError(f"unknown job type '{row[0]}'")
        except IndexError:
            print(f"Row {i} failed: missing columns")
            continue
        except ValueError as err:
            print(f"Row {i} failed: {err}")
            continue
        HISTORY.append(entry)
        print(f"Row {i}: {entry.operation} -> {entry.result}")


def generate_sample_csv(filename: str = "batch_sample.csv") -> None:
    rows = [
        ["# job_type", "args..."],
        ["mass_to_mole", "18.015", "18.015", "g", "mol"],
        ["mole_to_mass", "2.5", "32.0", "mol", "g"],
        ["ideal_gas", "V", "1.0", "atm", "0", "L", "2.0", "mol", "25", "C", "L"],
    ]
    with open(filename, "w", newline="", encoding="utf-8") as f:
        csv.writer(f).writerows(rows)
    print(f"Wrote {filename}")


def prompt_float(prompt: str) -> float:
    while True:
        try:
            return float(input(prompt))
        except ValueError:
            print("Enter a number.")


def prompt_quantity(label: str, units: str) -> Quantity:
    return prompt_float(f"{label} value: "), input(f"{label} unit ({units}): ").strip()


def prompt_target(prompt: str, choices: tuple[str, ...]) -> str:
    target = input(prompt).strip().upper()
    if target not in choices:
        raise ValueError(f"Choose one of: {', '.join(choices)}")
    return target


def prompt_molar_mass() -> float:
    while True:
        try:
            mass = calculate_molar_mass(input("Molar mass in g/mol (number or formula): "))
        except ValueError as err:
            print(err)
        else:
            print(f"Molar mass: {mass:.4f} g/mol")
            return mass


def prompt_environment(env_type: str) -> Quantity:
    label = env_type.capitalize()
    units = "C, K, F" if env_type == "temperature" else "atm, Pa, kPa, Torr"
    while True:
        raw = input(f"{label} value or preset (STP, RT, SATP): ").strip()
        try:
            if raw.upper() in PRESET_TEMPS_C:
                return parse_environment_input(raw, "", env_type)
            return float(raw), input(f"{label} unit ({units}): ").strip()
        except ValueError:
            print(f"Invalid {env_type}. Enter a number, STP, RT, or SATP.")


def main_menu() -> None:
    while True:
        print(f"\n{MENU}")
        choice = input("Select an option (1-12): ").strip()
        entry = None
        try:
            if choice == "1":
                sub = input("(1) Mass to moles, (2) Moles to mass: ").strip()
                if sub not in ("1", "2"):
                    raise ValueError("Choose 1 or 2.")
                value = prompt_float("Input value: ")
                src_unit = input("Input unit (e.g. g, kg, mol): ").strip()
                molar_mass = prompt_molar_mass()
                dst_unit = input("Output unit (e.g. mol, mmol, g): ").strip()
                mode = "mass_to_mole" if sub == "1" else "mole_to_mass"
                entry = calc_mass_mole_conversion(value, mode, molar_mass, src_unit, dst_unit)

            elif choice == "2":
                target = prompt_target("Solve for (P, V, n, T): ", ("P", "V", "N", "T"))
                pressure = UNKNOWN if target == "P" else prompt_environment("pressure")
                volume = UNKNOWN if target == "V" else prompt_quantity("Volume", "L, mL, m3")
                moles = UNKNOWN if target == "N" else prompt_quantity("Moles", "mol, mmol")
                temperature = UNKNOWN if target == "T" else prompt_environment("temperature")
                out_unit = input("Output unit: ").strip()
                entry = solve_ideal_gas(target, pressure, volume, moles, temperature, out_unit)

            elif choice == "3":
                mass = prompt_quantity("Mass", "g, kg, lb")
                molar_mass = prompt_molar_mass()
                pressure = prompt_environment("pressure")
                temperature = prompt_environment("temperature")
                out_unit = input("Output volume unit (L, mL, m3): ").strip()
                entry = calc_gas_volume_from_mass(mass, molar_mass, pressure, temperature, out_unit)

            elif choice == "4":
                category = input("Category ((M)ass, (V)olume, (P)ressure, (MOL)e): ")
                category = category.strip().upper()
                kind = "mole" if category.startswith("MOL") else {
                    "M": "mass", "V": "volume", "P": "pressure"}.get(category[:1])
                if kind is None:
                    raise ValueError("Unknown category.")
                value = prompt_float("Value: ")
                src = input("From unit: ").strip()
                dst = input("To unit: ").strip()
                print(f"{value:.10g} {src} = {convert(kind, value, src, dst):.10g} {dst}")

            elif choice == "5":
                process_batch_file(input("CSV file path: ").strip())

            elif choice == "6":
                generate_sample_csv()

            elif choice == "7":
                if not HISTORY:
                    print("No calculations yet.")
                for i, past in enumerate(HISTORY, start=1):
                    print(f"{i}. {past.operation} -> {past.result}")

            elif choice == "8":
                print(UNIT_HELP)

            elif choice == "9":
                return

            elif choice == "10":
                target = prompt_target("Solve for (M, n, V): ", ("M", "N", "V"))
                molarity = UNKNOWN if target == "M" else prompt_quantity("Molarity", "M, mol/L, mM")
                moles = UNKNOWN
                if target != "N":
                    source = input("Amount of substance: (1) moles, (2) mass + formula: ").strip()
                    if source == "1":
                        moles = prompt_quantity("Moles", "mol, mmol")
                    elif source == "2":
                        mass, mass_unit = prompt_quantity("Mass", "g, kg, ...")
                        molar_mass = prompt_molar_mass()
                        moles = (mass * unit_factor("mass", mass_unit) / molar_mass, "mol")
                    else:
                        raise ValueError("Choose 1 or 2.")
                volume = UNKNOWN if target == "V" else prompt_quantity("Volume", "L, mL, m3")
                out_unit = input("Output unit: ").strip()
                entry = calc_molarity(target, molarity, moles, volume, out_unit)

            elif choice == "11":
                target = prompt_target("Solve for (M1, V1, M2, V2): ", ("M1", "V1", "M2", "V2"))
                m1 = UNKNOWN if target == "M1" else prompt_quantity("M1", "M, mol/L, mM")
                v1 = UNKNOWN if target == "V1" else prompt_quantity("V1", "L, mL, m3")
                m2 = UNKNOWN if target == "M2" else prompt_quantity("M2", "M, mol/L, mM")
                v2 = UNKNOWN if target == "V2" else prompt_quantity("V2", "L, mL, m3")
                out_unit = input("Output unit: ").strip()
                entry = calc_dilution(target, m1, v1, m2, v2, out_unit)

            elif choice == "12":
                target = prompt_target("Solve for (d, m, V): ", ("D", "M", "V"))
                mass = UNKNOWN if target == "M" else prompt_quantity("Mass", "g, kg, ...")
                volume = UNKNOWN if target == "V" else prompt_quantity("Volume", "mL, L, m3")
                density = UNKNOWN if target == "D" else prompt_quantity("Density", "g/mL, kg/L")
                out_unit = input("Output unit: ").strip()
                entry = calc_liquid_density(target, mass, volume, density, out_unit)

            else:
                print("Invalid option. Choose 1-12.")

            if entry:
                HISTORY.append(entry)
                print(f"\n{entry}")
        except (ValueError, OSError) as err:
            print(f"Error: {err}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-b", "--batch", metavar="CSV", help="run the jobs in a CSV file and exit")
    args = parser.parse_args()

    try:
        if args.batch:
            process_batch_file(args.batch)
        else:
            main_menu()
    except (EOFError, KeyboardInterrupt):
        print()
    except (OSError, UnicodeDecodeError) as err:
        raise SystemExit(f"error: {err}")


if __name__ == "__main__":
    main()
