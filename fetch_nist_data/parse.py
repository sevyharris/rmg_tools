# Script to parse HTML data from NIST or other sources

from bs4 import BeautifulSoup
import pandas as pd
import re
from typing import List, Dict, Optional, Union


def parse_html_file(filename: str) -> BeautifulSoup:
    """
    Parse an HTML file and return a BeautifulSoup object.
    
    Args:
        filename: Path to the HTML file
        
    Returns:
        BeautifulSoup object for parsing
    """
    with open(filename, 'r', encoding='utf-8') as f:
        html_content = f.read()
    
    return BeautifulSoup(html_content, 'html.parser')


def extract_hf_data(soup: BeautifulSoup) -> pd.DataFrame:
    """
    Extract enthalpy of formation (Hf) data from NIST HTML page.
    
    Args:
        soup: BeautifulSoup object containing parsed HTML
        
    Returns:
        DataFrame with columns: Quantity, Value, Uncertainty, Units, Method, Reference, Comment
    """
    all_data = []
    
    # Find all tables
    for table in soup.find_all('table', class_='data'):
        rows = table.find_all('tr')
        
        for row in rows:
            cells = row.find_all(['td', 'th'])
            if len(cells) >= 3:
                # Check if this is a data row (not header)
                if cells[0].name == 'td':
                    quantity_text = cells[0].get_text(strip=True)
                    
                    # Look for enthalpy of formation entries
                    if 'ΔfH°' in quantity_text or 'fH°' in quantity_text or ('sub>f</sub>' in str(cells[0]) and 'H°' in str(cells[0])):
                        # Extract value with uncertainty
                        value_text = cells[1].get_text(strip=True)
                        
                        # Parse value and uncertainty
                        # Format can be: "-74.87" or "-74.6 ± 0.3"
                        value_match = re.match(r'([+-]?[0-9.]+)\s*(?:±\s*([0-9.]+))?', value_text)
                        if value_match:
                            value = float(value_match.group(1))
                            uncertainty = float(value_match.group(2)) if value_match.group(2) else None
                        else:
                            continue
                        
                        # Extract other fields
                        units = cells[2].get_text(strip=True) if len(cells) > 2 else ''
                        method = cells[3].get_text(strip=True) if len(cells) > 3 else ''
                        reference = cells[4].get_text(strip=True) if len(cells) > 4 else ''
                        comment = cells[5].get_text(strip=True) if len(cells) > 5 else ''
                        
                        all_data.append({
                            'Quantity': quantity_text,
                            'Value': value,
                            'Uncertainty': uncertainty,
                            'Units': units,
                            'Method': method,
                            'Reference': reference,
                            'Comment': comment
                        })
    
    # Create DataFrame
    if all_data:
        df = pd.DataFrame(all_data)
        return df
    else:
        return pd.DataFrame()


def extract_s_gas_data(soup: BeautifulSoup) -> pd.DataFrame:
    """
    Extract gas phase entropy (S°gas) data from NIST HTML page.
    
    Args:
        soup: BeautifulSoup object containing parsed HTML
        
    Returns:
        DataFrame with columns: Quantity, Value, Uncertainty, Units, Method, Reference, Comment
    """
    all_data = []
    
    # Find all tables
    for table in soup.find_all('table', class_='data'):
        rows = table.find_all('tr')
        
        for row in rows:
            cells = row.find_all(['td', 'th'])
            if len(cells) >= 3:
                # Check if this is a data row (not header)
                if cells[0].name == 'td':
                    quantity_text = cells[0].get_text(strip=True)
                    
                    # Look for entropy entries (S°gas)
                    if 'S°' in quantity_text and 'gas' in quantity_text:
                        # Extract value with uncertainty
                        value_text = cells[1].get_text(strip=True)
                        
                        # Parse value and uncertainty
                        # Format can be: "186.25" or "188.66 ± 0.42"
                        value_match = re.match(r'([+-]?[0-9.]+)\s*(?:±\s*([0-9.]+))?', value_text)
                        if value_match:
                            value = float(value_match.group(1))
                            uncertainty = float(value_match.group(2)) if value_match.group(2) else None
                        else:
                            continue
                        
                        # Extract other fields
                        units = cells[2].get_text(strip=True) if len(cells) > 2 else ''
                        method = cells[3].get_text(strip=True) if len(cells) > 3 else ''
                        reference = cells[4].get_text(strip=True) if len(cells) > 4 else ''
                        comment = cells[5].get_text(strip=True) if len(cells) > 5 else ''
                        
                        all_data.append({
                            'Quantity': quantity_text,
                            'Value': value,
                            'Uncertainty': uncertainty,
                            'Units': units,
                            'Method': method,
                            'Reference': reference,
                            'Comment': comment
                        })
    
    # Create DataFrame
    if all_data:
        df = pd.DataFrame(all_data)
        return df
    else:
        return pd.DataFrame()


def extract_cp_data(soup: BeautifulSoup) -> pd.DataFrame:
    """
    Extract heat capacity (Cp) data over temperature from NIST HTML page.
    
    Args:
        soup: BeautifulSoup object containing parsed HTML
        
    Returns:
        DataFrame with columns: Temperature (K), Cp (J/mol*K), Reference, Comment
    """
    all_data = []
    
    # Find all tables with Cp data
    # Look for h3 tags with "Constant pressure heat capacity"
    for h3 in soup.find_all('h3'):
        if 'Constant pressure heat capacity' in h3.get_text():
            # Find the table following this heading
            table = h3.find_next('table')
            if table:
                # Extract table rows
                rows = table.find_all('tr')
                
                # Get headers (first row)
                if len(rows) > 0:
                    # Skip header row and process data rows
                    for row in rows[1:]:
                        cells = row.find_all(['td', 'th'])
                        if len(cells) >= 2:
                            # Extract Cp value
                            cp_text = cells[0].get_text(strip=True)
                            # Extract temperature
                            temp_text = cells[1].get_text(strip=True)
                            
                            # Parse Cp value (remove uncertainty if present)
                            cp_match = re.match(r'([0-9.]+)', cp_text)
                            if cp_match:
                                cp_value = float(cp_match.group(1))
                            else:
                                continue
                            
                            # Parse temperature
                            temp_match = re.match(r'([0-9.]+)', temp_text)
                            if temp_match:
                                temp_value = float(temp_match.group(1))
                            else:
                                continue
                            
                            # Extract reference and comment if available
                            reference = cells[2].get_text(strip=True) if len(cells) > 2 else ''
                            comment = cells[3].get_text(strip=True) if len(cells) > 3 else ''
                            
                            all_data.append({
                                'Temperature (K)': temp_value,
                                'Cp (J/mol*K)': cp_value,
                                'Reference': reference,
                                'Comment': comment
                            })
    
    # Create DataFrame
    if all_data:
        df = pd.DataFrame(all_data)
        # Sort by temperature
        df = df.sort_values('Temperature (K)').reset_index(drop=True)
        return df
    else:
        return pd.DataFrame()


def extract_all_tables(soup: BeautifulSoup) -> List[pd.DataFrame]:
    """
    Extract all HTML tables from a BeautifulSoup object into pandas DataFrames.
    
    Args:
        soup: BeautifulSoup object containing parsed HTML
        
    Returns:
        List of pandas DataFrames, one for each table found
    """
    tables = soup.find_all('table')
    dataframes = []
    
    for table in tables:
        # Try to extract headers
        headers = []
        header_row = table.find('tr')
        if header_row:
            headers = [th.get_text(strip=True) for th in header_row.find_all(['th', 'td'])]
        
        # Extract all rows
        rows = []
        for tr in table.find_all('tr')[1:] if headers else table.find_all('tr'):
            cells = [td.get_text(strip=True) for td in tr.find_all(['td', 'th'])]
            if cells:
                rows.append(cells)
        
        # Create DataFrame
        if rows:
            if headers and len(headers) == len(rows[0]):
                df = pd.DataFrame(rows, columns=headers)
            else:
                df = pd.DataFrame(rows)
            dataframes.append(df)
    
    return dataframes


def extract_metadata(soup: BeautifulSoup) -> Dict[str, str]:
    """
    Extract metadata from the HTML page (title, formula, molecular weight, etc.).
    
    Args:
        soup: BeautifulSoup object containing parsed HTML
        
    Returns:
        Dictionary of metadata key-value pairs
    """
    metadata = {}
    
    # Extract title
    title = soup.find('h1', id='Top')
    if title:
        metadata['title'] = title.get_text(strip=True)
    
    # Extract list items that contain metadata
    for li in soup.find_all('li'):
        # Look for strong tags that often contain labels
        strong = li.find('strong')
        if strong:
            label = strong.get_text(strip=True).rstrip(':')
            # Get the rest of the text content
            text = li.get_text(strip=True)
            # Remove the label from the text
            value = text.replace(label + ':', '', 1).strip()
            if value and not value.startswith('http'):  # Skip links
                metadata[label] = value
    
    return metadata


def parse_nist_page(filename: str = 'output.html') -> Dict[str, Union[pd.DataFrame, Dict, List]]:
    """
    Main function to parse a NIST HTML page.
    
    Args:
        filename: Path to the HTML file
        
    Returns:
        Dictionary containing:
            - 'metadata': Dictionary of metadata
            - 'hf_data': DataFrame with enthalpy of formation data
            - 's_gas_data': DataFrame with gas phase entropy data
            - 'cp_data': DataFrame with Cp vs Temperature data
            - 'all_tables': List of all DataFrames from all tables
    """
    print(f"Parsing {filename}...")
    
    # Parse HTML
    soup = parse_html_file(filename)
    
    # Extract metadata
    metadata = extract_metadata(soup)
    print(f"\nFound metadata fields: {len(metadata)}")
    if 'title' in metadata:
        print(f"Compound: {metadata['title']}")
    
    # Extract Hf data specifically
    hf_df = extract_hf_data(soup)
    print(f"Extracted {len(hf_df)} Hf data points")
    
    # Extract S gas data specifically
    s_gas_df = extract_s_gas_data(soup)
    print(f"Extracted {len(s_gas_df)} S gas data points")
    
    # Extract Cp data specifically
    cp_df = extract_cp_data(soup)
    print(f"Extracted {len(cp_df)} Cp data points")
    
    # Extract all tables
    all_tables = extract_all_tables(soup)
    print(f"Found {len(all_tables)} total table(s)")
    
    result = {
        'metadata': metadata,
        'hf_data': hf_df,
        's_gas_data': s_gas_df,
        'cp_data': cp_df,
        'all_tables': all_tables
    }
    
    return result


def save_to_csv(df: pd.DataFrame, filename: str) -> None:
    """
    Save a DataFrame to a CSV file.
    
    Args:
        df: DataFrame to save
        filename: Output CSV filename
    """
    if df is not None and not df.empty:
        df.to_csv(filename, index=False)
        print(f"\nData saved to {filename}")
    else:
        print("No data to save")


if __name__ == "__main__":
    import sys
    
    # Get filename from command line or use default
    filename = sys.argv[1] if len(sys.argv) > 1 else 'output.html'
    
    # Parse the HTML file
    result = parse_nist_page(filename)
    
    # Display metadata
    print("\n" + "="*60)
    print("METADATA:")
    print("="*60)
    for key, value in result['metadata'].items():
        display_value = str(value)[:100] + "..." if len(str(value)) > 100 else str(value)
        print(f"{key}: {display_value}")
    
    # Display Hf data
    if not result['hf_data'].empty:
        print("\n" + "="*60)
        print("ENTHALPY OF FORMATION (Hf) DATA:")
        print("="*60)
        print(result['hf_data'].to_string(index=False))
        
        # Save to CSV
        save_to_csv(result['hf_data'], 'hf_data.csv')
        
        # Display statistics
        print("\n" + "="*60)
        print("Hf STATISTICS:")
        print("="*60)
        print(f"Number of entries: {len(result['hf_data'])}")
        print(f"Value range: {result['hf_data']['Value'].min():.2f} - {result['hf_data']['Value'].max():.2f} {result['hf_data']['Units'].iloc[0] if len(result['hf_data']) > 0 else ''}")
    else:
        print("\nNo Hf data found in the HTML file.")
    
    # Display S gas data
    if not result['s_gas_data'].empty:
        print("\n" + "="*60)
        print("GAS PHASE ENTROPY (S°gas) DATA:")
        print("="*60)
        print(result['s_gas_data'].to_string(index=False))
        
        # Save to CSV
        save_to_csv(result['s_gas_data'], 's_gas_data.csv')
        
        # Display statistics
        print("\n" + "="*60)
        print("S°gas STATISTICS:")
        print("="*60)
        print(f"Number of entries: {len(result['s_gas_data'])}")
        print(f"Value range: {result['s_gas_data']['Value'].min():.2f} - {result['s_gas_data']['Value'].max():.2f} {result['s_gas_data']['Units'].iloc[0] if len(result['s_gas_data']) > 0 else ''}")
    else:
        print("\nNo S gas data found in the HTML file.")
    
    # Display Cp data
    if not result['cp_data'].empty:
        print("\n" + "="*60)
        print("HEAT CAPACITY (Cp) DATA:")
        print("="*60)
        print(result['cp_data'].to_string(index=False))
        
        # Save to CSV
        save_to_csv(result['cp_data'], 'cp_data.csv')
        
        # Display statistics
        print("\n" + "="*60)
        print("Cp STATISTICS:")
        print("="*60)
        print(f"Temperature range: {result['cp_data']['Temperature (K)'].min():.1f} - {result['cp_data']['Temperature (K)'].max():.1f} K")
        print(f"Cp range: {result['cp_data']['Cp (J/mol*K)'].min():.2f} - {result['cp_data']['Cp (J/mol*K)'].max():.2f} J/mol*K")
    else:
        print("\nNo Cp data found in the HTML file.")
