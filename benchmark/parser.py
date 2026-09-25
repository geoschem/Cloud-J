"""J-value output parser for Cloud-J benchmark harness.

Parses the "Fast-J ----J-values----" output blocks from both the
Fortran and C++ Cloud-J standalone executables into structured data.
"""

from typing import List

from benchmark.types import JValueBlock


_MARKER = "Fast-J ----J-values----"


def parse_jvalues(stdout: str) -> List[JValueBlock]:
    """Parse all J-value blocks from executable stdout.

    Looks for 'Fast-J ----J-values----' markers, then reads the species header
    line (starting with 'L=') and subsequent numeric rows.

    Each row has format: layer_num val1 val2 val3 ...
    where values are in Fortran e9.2 scientific notation (e.g., '2.54E-09').

    Args:
        stdout: Full captured stdout from a Cloud-J executable run.

    Returns:
        List of JValueBlock objects, one per SZA block found in the output.
    """
    blocks: List[JValueBlock] = []
    lines = stdout.splitlines()
    num_lines = len(lines)
    i = 0
    sza_index = 0

    while i < num_lines:
        # Look for the marker line
        if _MARKER in lines[i]:
            # Next line should be the species header starting with "L="
            i += 1
            if i >= num_lines:
                break

            header_line = lines[i].strip()
            if not header_line.startswith("L="):
                # Unexpected format; skip this marker
                continue

            # Parse species names from header
            # Format: "L=  O2       O3       O3(1D)   NO   ..."
            # or:     "L=      O2       O3   O3(1D)"
            header_content = header_line[2:]  # Strip "L="
            species = header_content.split()

            # Parse numeric rows
            values = {}
            i += 1
            while i < num_lines:
                row_line = lines[i].strip()
                if not row_line:
                    # Empty line ends the block
                    break

                parts = row_line.split()
                if not parts:
                    break

                # First element should be the layer number (integer)
                try:
                    layer_num = int(parts[0])
                except ValueError:
                    # Not a numeric row; block has ended
                    break

                # Remaining elements are J-values in scientific notation
                row_values: List[float] = []
                for val_str in parts[1:]:
                    try:
                        row_values.append(float(val_str))
                    except ValueError:
                        # Non-numeric value encountered; stop parsing this row
                        break

                if row_values:
                    values[layer_num] = row_values

                i += 1

            block = JValueBlock(
                sza_index=sza_index,
                species=species,
                values=values,
            )
            blocks.append(block)
            sza_index += 1
        else:
            i += 1

    return blocks
