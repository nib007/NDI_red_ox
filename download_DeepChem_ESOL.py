
import pandas as pd
import requests
from io import StringIO

# Alternative ESOL dataset URL (Zenodo mirror)
url = "https://zenodo.org/record/3778405/files/esol.csv?download=1"

# Download the file
response = requests.get(url)

if response.status_code == 200:
    # Convert response text to a DataFrame
    esol_df = pd.read_csv(StringIO(response.text))
    print("Download successful!")
    print(esol_df.head())  # Show first 5 rows
else:
    print(f"Failed to download ESOL dataset. Status code: {response.status_code}")
